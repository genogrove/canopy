# SPDX-License-Identifier: GPL-3.0-or-later
"""The `text` table is for a human to read; `tsv`/`json`/`bed` are for a parser.

The two are allowed to diverge, and these tests pin *where*: text drops columns that carry no
per-row information and abbreviates nested cells, while the machine formats stay a stable
rectangle with every field intact.
"""

from __future__ import annotations

import json

from genogrove_canopy.cli import _render

# An enhancer answer, trimmed to three rows. Five of its columns are identical on every row.
_ROWS = [
    {"chrom": "chrX", "start": 67544302, "end": 67545353, "strand": ".",
     "name": "enh:promoter->AR", "type": "enhancer", "class": "promoter", "score": 0.99999,
     "n": 1, "cohort": "EFO:0005726", "target": "AR",
     "ccre_overlap": [{"id": "EH38E3936120", "class": "pELS", "bp": 332},
                      {"id": "EH38E3936121", "class": "PLS", "bp": 277}]},
    {"chrom": "chrX", "start": 67551264, "end": 67551763, "strand": ".",
     "name": "enh:genic->AR", "type": "enhancer", "class": "genic", "score": 0.6421,
     "n": 1, "cohort": "EFO:0005726", "target": "AR", "ccre_overlap": []},
    {"chrom": "chrX", "start": 67538334, "end": 67538833, "strand": ".",
     "name": "enh:intergenic->AR", "type": "enhancer", "class": "intergenic", "score": 0.47963,
     "n": 1, "cohort": "EFO:0005726", "target": "AR",
     "ccre_overlap": [{"id": "EH38E4521690", "class": "TF", "bp": 168}]},
]
_STDOUT = "\n".join(json.dumps(r) for r in _ROWS)


def test_constant_columns_are_stated_once_not_repeated() -> None:
    """`cohort`, `target`, `type`, `strand` and `n` are identical on all three rows."""
    lines = _render(_STDOUT, "text").strip().splitlines()
    shared, header = lines[0], lines[1]
    for col in ("cohort=EFO:0005726", "target=AR", "type=enhancer", "n=1", "strand=."):
        assert col in shared
        assert col.split("=")[0] not in header.split(), f"{col} was collapsed but still a column"
    # the columns that actually differ must survive
    for col in ("start", "end", "class", "score", "ccre_overlap"):
        assert col in header.split()


def test_nested_cells_are_abbreviated_not_repr_dumped() -> None:
    text = _render(_STDOUT, "text")
    assert "pELS:332 PLS:277" in text
    assert "EH38E3936120" not in text, "ids belong in --format json, not the human table"
    assert "{'id'" not in text and '{"id"' not in text


def test_empty_list_reads_as_a_finding_not_a_blank() -> None:
    """`ccre_overlap: []` means *no overlap*, which is a result. A blank cell reads as missing."""
    row = [ln for ln in _render(_STDOUT, "text").splitlines() if "genic" in ln and "inter" not in ln]
    assert row and row[0].rstrip().endswith("-")


def test_machine_formats_keep_every_column_and_the_full_structure() -> None:
    """A parser cannot know a column vanished because it happened to be constant."""
    header = _render(_STDOUT, "tsv").splitlines()[0].split("\t")
    for col in ("chrom", "strand", "type", "n", "cohort", "target", "ccre_overlap"):
        assert col in header

    first = json.loads(_render(_STDOUT, "json").splitlines()[0])
    assert first["ccre_overlap"][0]["id"] == "EH38E3936120", "json must keep the cCRE ids"
    assert first["cohort"] == "EFO:0005726"


def test_a_single_record_table_is_unchanged() -> None:
    """Collapsing needs something to compare against; one row has no constant columns by
    definition, and hoisting all of them would leave an empty table."""
    one = json.dumps(_ROWS[0])
    lines = _render(one, "text").strip().splitlines()
    assert lines[0].split()[:3] == ["chrom", "start", "end"], "no shared-values line for one row"
