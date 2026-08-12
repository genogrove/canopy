# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure codegen helpers — no API key or pygenogrove needed."""

import json

from canopy import llm
from canopy.cli import _render


def test_build_system_prompt_injects_resources_block():
    prompt = llm.build_system_prompt("- `GENCODE_HUMAN` (str): a grove path")
    assert "## Available resources" in prompt
    assert "GENCODE_HUMAN" in prompt
    assert "## The GENCODE Grove model" in prompt  # earlier sections preserved
    assert "TODO: injected at runtime" not in prompt  # placeholder dropped
    assert prompt.index("GENCODE Grove model") < prompt.index("Available resources")


_REC = '{"chrom": "chr7", "start": 100, "end": 200, "name": "EGFR", "strand": "+"}'


def test_render_bed_converts_to_half_open():
    out = _render(_REC, "bed")
    assert out.startswith("#chrom\tstart\tend\tname\tscore\tstrand\n")
    assert "chr7\t100\t201\tEGFR\t.\t+" in out  # end 200 -> 201, default score "."


def test_render_tsv_and_json():
    assert "chrom\tstart\tend\tname\tstrand" in _render(_REC, "tsv")
    assert json.loads(_render(_REC, "json").strip())["name"] == "EGFR"  # grove-native, unconverted


def test_render_scalar_passes_through():
    assert _render("count: 42", "bed").strip() == "count: 42"


def test_strip_code_fence():
    assert llm._strip_code_fence("```python\nprint(1)\n```") == "print(1)\n"
    assert llm._strip_code_fence("```\nx = 1\n```") == "x = 1\n"
    assert llm._strip_code_fence("print(2)") == "print(2)\n"  # unfenced passes through
