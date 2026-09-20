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


def test_bare_cohort_line_does_not_capture_the_next_line():
    """A question naming no biosample yields a bare ``COHORT:`` — which must parse as no cohort.

    Regression: the pattern used ``\\s*``, which matches newlines, so it stepped over the empty
    value and captured the ``TARGETS:`` line as the cohort name. The host then resolved that
    against the ENCODE catalog. Found by benchmarks/bench.py on both models at once — the
    give-away that it was canopy's parsing, not the model's declaration.
    """
    cohort, targets, _ = llm.parse_targets_and_code(
        "COHORT:\nTARGETS: []\n\n```python\nprint(1)\n```"
    )
    assert cohort == ""
    assert targets == []

    # The array on its own line is the same trap one line further down.
    cohort, _, _ = llm.parse_targets_and_code(
        "COHORT:\nTARGETS:\n[]\n\n```python\nprint(1)\n```"
    )
    assert cohort == ""

    # ...while a real cohort still parses, quoted or bare.
    for line in ("COHORT: MCF-7", "COHORT: 'MCF-7'"):
        cohort, _, _ = llm.parse_targets_and_code(f"{line}\nTARGETS: []\n\n```python\nprint(1)\n```")
        assert cohort == "MCF-7", line


def test_placeholder_cohorts_normalise_to_no_cohort():
    """`COHORT: none` means no cohort — not a biosample literally called "none".

    system.md asks the model to omit the line entirely, but models routinely write a placeholder
    instead, and the host would hand that to the ENCODE catalog. Normalising here rather than
    hardening the prompt keeps canopy model-agnostic: scoring one vendor's phrasing as the only
    correct one measures a formatting habit, not capability.
    """
    for value in ("none", "None", "(none)", "'none'", "N/A", "null", "-", "not specified"):
        cohort, _, _ = llm.parse_targets_and_code(
            f"COHORT: {value}\nTARGETS: []\n\n```python\nprint(1)\n```"
        )
        assert cohort == "", f"{value!r} should mean no cohort, got {cohort!r}"

    # A real biosample whose name merely contains a sentinel substring is untouched.
    for value in ("NALM-6", "HL-60", "MCF-7"):
        cohort, _, _ = llm.parse_targets_and_code(
            f"COHORT: {value}\nTARGETS: []\n\n```python\nprint(1)\n```"
        )
        assert cohort == value, value


def test_strip_code_fence():
    assert llm._strip_code_fence("```python\nprint(1)\n```") == "print(1)\n"
    assert llm._strip_code_fence("```\nx = 1\n```") == "x = 1\n"
    assert llm._strip_code_fence("print(2)") == "print(2)\n"  # unfenced passes through
