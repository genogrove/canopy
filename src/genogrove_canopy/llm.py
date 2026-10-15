# SPDX-License-Identifier: GPL-3.0-or-later
"""The agent loop, via the Anthropic API.

Given a natural-language question, the model drives a ``pygenogrove`` grove through the tools
the host provides — ``attach`` (put a pinned cohort's layers into the grove), ``build`` (run a
program that inserts its own data into the grove) and ``run`` (execute a query program against
it) — reads what comes back, and calls again or answers. The system prompt in
``prompts/system.md`` is its API contract; the resource context (see
:mod:`genogrove_canopy.resources`) is injected into it.

The generated code is untrusted: it must only ever run through :mod:`genogrove_canopy.sandbox`,
never ``exec``'d here. The host-side tool handlers (``cli.Session``) do exactly that.
"""

from __future__ import annotations

import re
from pathlib import Path

DEFAULT_MODEL = "claude-opus-4-8"

#: Hard cap on model turns per question. The model gets this many replies; a reply with no tool
#: call is the answer. Running out means abstain — the trace says what was tried.
MAX_TURNS = 8

_SYSTEM_MD = Path(__file__).with_name("prompts") / "system.md"

#: The per-question layers the host knows how to attach (see `cli.Session.attach`).
_LAYERS = ("enhancers", "sv", "hichip")

#: Ways a model spells "no biosample" instead of leaving the cohort empty. Normalised to "" so
#: the host doesn't hand a placeholder to the ENCODE catalog. Deliberately permissive: canopy is
#: model-agnostic, and treating one vendor's phrasing as the only valid one would score every
#: other model down for a formatting habit, not a capability.
_NO_COHORT = frozenset((
    "", "none", "n/a", "na", "null", "nil", "-", "--",
    "unspecified", "not specified", "not applicable",
))
# The system prompt's runtime-injected datasets block (the TODO placeholder).
_RESOURCES_HEADING = "## Available resources"

#: The agent's action space: what the model can do to the grove, and nothing else.
TOOLS = [
    {
        "name": "attach",
        "description": (
            "Attach a cohort's per-question layers into GROVE before querying them: rE2G "
            "enhancer->gene links, PCAWG structural variants and/or TCGA HiChIP contacts. "
            "Returns what is now attached and the COHORTS / SV_COHORTS / HICHIP_COHORTS your "
            "code must filter on. Call it once per cohort set; the grove keeps what it holds."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "layers": {
                    "type": "array",
                    "items": {"type": "string", "enum": list(_LAYERS)},
                    "description": "Exactly the layers the question needs.",
                },
                "cohort": {
                    "type": "string",
                    "description": (
                        "The tissue / disease / cell line the question implies, as one of the "
                        "terms under 'Available resources'; several `;`-separated for a "
                        "comparison. Empty if the question names none."
                    ),
                },
            },
            "required": ["layers"],
        },
    },
    {
        "name": "build",
        "description": (
            "Run one program with GROVE **mutable** (insert, add_edge, add_external_key) to put "
            "your own data into the grove — the user's files listed under 'Available resources', "
            "read with open(path). What you insert stays for the session, like an attached "
            "cohort; the next `run` queries it. Returns the program's stdout (or the error). Not "
            "for querying: print only what you need to confirm the build."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"code": {"type": "string", "description": "The Python program."}},
            "required": ["code"],
        },
    },
    {
        "name": "run",
        "description": (
            "Run one self-contained pygenogrove program against GROVE (read-only) in the sandbox "
            "and return its stdout (or the error). Print records as the Rules say; read them; "
            "run again or answer."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"code": {"type": "string", "description": "The Python program."}},
            "required": ["code"],
        },
    },
]


def build_system_prompt(resources_block: str) -> str:
    """The system prompt: ``system.md`` with the resources block injected.

    ``resources_block`` replaces the placeholder under "Available resources" — it
    names the variables holding each dataset path and what they are. The generated
    code emits canonical records (see the Rules); the host renders the user's
    ``--format`` choice, so the format is not part of this prompt.
    """
    text = _SYSTEM_MD.read_text(encoding="utf-8")
    head, _, _tail = text.partition(_RESOURCES_HEADING)
    return f"{head}{_RESOURCES_HEADING}\n\n{resources_block.strip()}\n"


def answer(question: str, system_prompt: str, tools, *, model: str = DEFAULT_MODEL,
           max_turns: int = MAX_TURNS, on_turn=None, on_usage=None):
    """Answer ``question`` by letting the model drive ``tools`` until it replies without one.

    ``tools`` is the host side of :data:`TOOLS`: an object with ``attach(layers, cohort="")``,
    ``build(code)`` and ``run(code)``, each returning the text the model reads next
    (``cli.Session``). Nothing is executed here.

    Returns ``(text, trace)``: ``text`` is the model's final reply, ``""`` if it ran out of
    turns (abstain); ``trace`` is every tool call in order as ``{"tool", "input", "output"}`` —
    the provenance of the answer. ``on_turn(entry)`` is called as each is made. ``on_usage``
    receives each response's ``usage`` (cost accounting for a benchmark harness, not the
    package). Raises ``RuntimeError`` if the model declines.
    """
    import anthropic  # lazy: keeps the module importable without the SDK/key

    client = anthropic.Anthropic()  # reads ANTHROPIC_API_KEY from the environment
    messages = [{"role": "user", "content": question}]
    trace = []
    for _ in range(max_turns):
        response = client.messages.create(
            model=model,
            max_tokens=16000,  # room for adaptive thinking + a program; non-streaming-safe
            thinking={"type": "adaptive"},
            system=system_prompt,
            tools=TOOLS,
            messages=messages,
        )
        if on_usage is not None:
            on_usage(response.usage)
        if response.stop_reason == "refusal":
            raise RuntimeError("the model declined to answer this question")
        messages.append({"role": "assistant", "content": response.content})
        calls = [b for b in response.content if b.type == "tool_use"]
        if not calls:
            return "".join(b.text for b in response.content if b.type == "text").strip(), trace
        results = []
        for call in calls:
            output = _dispatch(tools, call.name, call.input)
            entry = {"tool": call.name, "input": call.input, "output": output}
            trace.append(entry)
            if on_turn is not None:
                on_turn(entry)
            results.append({"type": "tool_result", "tool_use_id": call.id, "content": output})
        messages.append({"role": "user", "content": results})
    return "", trace


def _dispatch(tools, name: str, args: dict) -> str:
    """Route one tool call to the host, normalising what a model may send: a fenced program, a
    placeholder cohort, a layer the host doesn't have. A bad call is reported back to the model
    as text, never raised — the model is the one that can fix it."""
    if name in ("run", "build"):
        return getattr(tools, name)(_strip_code_fence(str(args.get("code", ""))))
    if name == "attach":
        cohort = str(args.get("cohort") or "").strip().strip("\"'")
        if cohort.lower() in _NO_COHORT:
            cohort = ""
        layers = [w for w in (str(x).strip().lower() for x in args.get("layers") or []) if w in _LAYERS]
        return tools.attach(layers, cohort)
    return f"error: unknown tool {name!r}"


def _strip_code_fence(text: str) -> str:
    """Return the Python program from a reply or tool input.

    Prefer an explicit ```python fence; else the last ``` block; else the text as-is. The tool
    takes bare code, but be robust to a model that fences it anyway.
    """
    m = re.search(r"```python\n(.*?)```", text, re.DOTALL)
    if not m:
        blocks = re.findall(r"```(?:\w+)?\n(.*?)```", text, re.DOTALL)
        return (blocks[-1] if blocks else text).strip() + "\n"
    return m.group(1).strip() + "\n"
