# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure codegen helpers — no API key or pygenogrove needed."""

import json

from genogrove_canopy import llm
from genogrove_canopy.cli import _render


def test_build_system_prompt_injects_resources_block():
    prompt = llm.build_system_prompt("- `GENCODE_HUMAN` (str): a grove path")
    assert "## Available resources" in prompt
    assert "GENCODE_HUMAN" in prompt
    assert "## The GENCODE Grove model" in prompt  # earlier sections preserved
    assert "TODO: injected at runtime" not in prompt  # placeholder dropped
    assert prompt.index("## The GENCODE Grove model") < prompt.index("## Available resources")


_REC = '{"chrom": "chr7", "start": 100, "end": 200, "name": "EGFR", "strand": "+"}'


def test_render_bed_converts_to_half_open():
    out = _render(_REC, "bed")
    assert out.startswith("#chrom\tstart\tend\tname\tscore\tstrand\n")
    assert "chr7\t100\t201\tEGFR\t.\t+" in out  # end 200 -> 201, default score "."


def test_render_tsv_and_json():
    # column order is deliberate (_record_columns): coordinates, then identity, then evidence
    assert "chrom\tstart\tend\tstrand\tname" in _render(_REC, "tsv")
    assert json.loads(_render(_REC, "json").strip())["name"] == "EGFR"  # grove-native, unconverted


def test_render_null_field_as_dot_not_python_none():
    """A null (e.g. a translocation's `size`) is `.` in tsv and text, never the literal `None`;
    json keeps `null` (#37)."""
    rec = '{"chrom": "chr8", "start": 1, "end": 1, "name": "x", "size": null}'
    for fmt in ("tsv", "text"):
        out = _render(rec, fmt)
        assert "None" not in out and out.rstrip().endswith(".")
    assert json.loads(_render(rec, "json").strip())["size"] is None
    # bed: an explicit null score/strand is "." like an absent one, never a malformed "None"
    assert _render(rec.replace('"size": null', '"score": null, "strand": null'), "bed").strip().endswith("x\t.\t.")
    # text: a column that is null in every row collapses to `size=.` above the table
    assert "size=." in _render("\n".join([rec] * 3), "text")


def test_render_scalar_passes_through():
    assert _render("count: 42", "bed").strip() == "count: 42"


def test_strip_code_fence():
    assert llm._strip_code_fence("```python\nprint(1)\n```") == "print(1)\n"
    assert llm._strip_code_fence("```\nx = 1\n```") == "x = 1\n"
    assert llm._strip_code_fence("print(2)") == "print(2)\n"  # unfenced passes through


class _Tools:
    """A host that records what the model asked for and answers canned text."""

    def __init__(self):
        self.calls = []

    def attach(self, layers, cohort=""):
        self.calls.append(("attach", layers, cohort))
        return "attached"

    def run(self, code):
        self.calls.append(("run", code))
        return "count: 1\n" if "ok" in code else "error:\nNameError: x"

    def build(self, code):
        self.calls.append(("build", code))
        return "built"


def test_dispatch_normalises_what_a_model_sends():
    """`attach` placeholders and `run` fences never reach the host as-is — "none" is no cohort
    (a biosample literally called "none" doesn't exist; the catalog would be searched for it),
    an unknown layer is dropped, a fenced program is unwrapped. Model-agnostic on purpose:
    scoring one vendor's formatting habit as the only correct one measures nothing."""
    t = _Tools()
    llm._dispatch(t, "attach", {"layers": ["Enhancers", "bogus"], "cohort": "'none'"})
    llm._dispatch(t, "run", {"code": "```python\nprint(1)\n```"})
    llm._dispatch(t, "build", {"code": "```\nins()\n```"})
    assert t.calls == [("attach", ["enhancers"], ""), ("run", "print(1)\n"), ("build", "ins()\n")]
    assert [tool["name"] for tool in llm.TOOLS] == ["attach", "build", "run"]
    for value in ("none", "N/A", "null", "-", "not specified"):
        assert llm._dispatch(_Tools(), "attach", {"layers": ["sv"], "cohort": value}) == "attached"
    assert llm._dispatch(t, "attach", {"layers": ["sv"], "cohort": "MCF-7"}) and t.calls[-1][2] == "MCF-7"
    assert llm._dispatch(t, "frobnicate", {}).startswith("error:")


def _scripted_client(monkeypatch, turns):
    """Fake `anthropic.Anthropic` whose `messages.create` returns one scripted reply per call.
    Each turn is a list of content blocks: ("text", s) or ("tool_use", name, input)."""
    import anthropic
    from types import SimpleNamespace as NS

    seen = []

    def create(**kw):
        seen.append(kw)
        blocks = []
        for i, b in enumerate(turns.pop(0)):
            if b[0] == "text":
                blocks.append(NS(type="text", text=b[1]))
            else:
                blocks.append(NS(type="tool_use", id=f"t{len(seen)}_{i}", name=b[1], input=b[2]))
        return NS(content=blocks, stop_reason="end_turn", usage=None)

    monkeypatch.setattr(anthropic, "Anthropic", lambda: NS(messages=NS(create=create)))
    return seen


def test_answer_loops_until_the_model_replies_without_a_tool(monkeypatch):
    """The agentic contract: a tool call's result goes back to the model as a `tool_result`,
    the model may call again (here: fix a failed program), and the reply with no tool call is
    the answer. The trace holds every call in order — the provenance."""
    seen = _scripted_client(monkeypatch, [
        [("tool_use", "attach", {"layers": ["sv"], "cohort": "breast"}),
         ("tool_use", "run", {"code": "print(x)\n"})],
        [("text", "that failed, retrying"), ("tool_use", "run", {"code": "ok()\n"})],
        [("text", "One rearrangement.")],
    ])
    tools = _Tools()
    text, trace = llm.answer("q", "sys", tools, model="m")
    assert text == "One rearrangement."
    assert [t["tool"] for t in trace] == ["attach", "run", "run"]
    assert trace[1]["output"].startswith("error:") and trace[2]["output"] == "count: 1\n"
    # the second request carries the first turn's results back, one per call, by id
    results = seen[1]["messages"][2]["content"]
    assert [r["type"] for r in results] == ["tool_result", "tool_result"]
    assert results[1]["content"].startswith("error:")
    assert seen[0]["tools"] is llm.TOOLS and seen[0]["system"] == "sys"


def test_answer_abstains_at_the_turn_cap(monkeypatch):
    """A model that never stops calling tools gets cut off: empty text, the trace intact."""
    _scripted_client(monkeypatch, [[("tool_use", "run", {"code": "ok()\n"})]] * 3)
    text, trace = llm.answer("q", "sys", _Tools(), model="m", max_turns=3)
    assert text == "" and len(trace) == 3
