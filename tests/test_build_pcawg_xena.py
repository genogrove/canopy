# SPDX-License-Identifier: GPL-3.0-or-later
"""The build tool's two judgment calls, on synthetic inputs: which specimen a tumour's RNA-seq
or copy number is joined to (the SV layer's WGS aliquot), and when a lifted segment is kept."""

import importlib.util
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "build_pcawg_xena", Path(__file__).parents[1] / "tools" / "build_pcawg_xena.py")
build = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build)


def _row(donor, specimen, aliquot, kind="Primary tumour - solid tissue", strategy="WGS"):
    return {"icgc_donor_id": donor, "icgc_specimen_id": specimen, "aliquot_id": aliquot,
            "dcc_specimen_type": kind, "library_strategy": strategy}


def test_resolve_specimens_prefers_the_specimens_own_sv_aliquot_then_the_donors_only_one():
    sv = {"W1": "LIRI-JP", "W2": "LIRI-JP", "W3": "LIRI-JP"}
    rows = [
        _row("D1", "SP1", "W1"),                            # tumour WGS in the SV set
        _row("D1", "SP1", "R1", strategy="RNA-Seq"),        # RNA under the same specimen: own aliquot
        _row("D1", "SPN", "N1", kind="Normal - blood derived"),
        _row("D2", "SP2", "W2"),
        _row("D2", "SP2b", "R2", strategy="RNA-Seq"),       # RNA under another specimen id, one SV tumour: donor fallback
        _row("D2", "SPN2", "N2", kind="Normal - blood derived"),
        _row("D2", "SPN2", "RN2", kind="Normal - blood derived", strategy="RNA-Seq"),  # normal RNA: never joined
        _row("D3", "SP3", "W3"),
        _row("D3", "SP3x", "W9"),                           # a second tumour WGS not in the SV set
        _row("D3", "SP3b", "R3", strategy="RNA-Seq"),       # one SV aliquot for the donor: still unambiguous
        _row("D4", "SP4", "W4"),                            # donor with no SV aliquot at all
        _row("D4", "SP4b", "R4", strategy="RNA-Seq"),
    ]
    out, own = build.resolve_specimens(rows, sv)
    assert out == {"SP1": "W1", "SP2": "W2", "SP2b": "W2", "SP3": "W3", "SP3b": "W3", "SP3x": "W3"}
    assert own == {"SP1", "SP2", "SP3"}
    # a donor with two SV tumours gets no guess for an unmatched specimen
    sv2 = {**sv, "W9": "LIRI-JP"}
    out2, _ = build.resolve_specimens(rows, sv2)
    assert "SP3b" not in out2 and out2["SP3x"] == "W9"


def test_pick_columns_keeps_one_column_per_aliquot_preferring_the_specimens_own():
    """OV-AU has donors with a primary (own specimen) and recurrences (other specimens): the
    fallback must not add a second column for the same SV tumour, and a donor with several
    fallback candidates and no own column gets none — the sequenced tumour is unknown."""
    resolve = {"SP1": "W1", "SP1r": "W1", "SP2b": "W2", "SP3m1": "W3", "SP3m2": "W3"}
    header = ["feature", "SP1r", "SP1", "SP2b", "SP3m1", "SP3m2", "SPunknown"]
    chosen, unmatched, dropped = build.pick_columns(header, resolve, own={"SP1"})
    assert chosen == {"W1": 2, "W2": 3}          # own column wins over the recurrence; single fallback kept
    assert unmatched == 1 and dropped == 3        # SP1r duplicate; SP3m1/SP3m2 ambiguous


class _Chain:
    """A fake liftover: shifts positions by +10 on chr1, drops a hole, flips a block on chr2,
    stretches chr3 by inserting 5 Mb after position 1,000,000."""

    def convert_coordinate(self, chrom, pos):
        if chrom == "chr1":
            if 500 <= pos < 600:
                return []
            return [(chrom, pos + 10, "+")]
        if chrom == "chr2":
            return [(chrom, pos, "-")]
        if chrom == "chr3":
            return [(chrom, pos + (5_000_000 if pos >= 1_000_000 else 0), "+")]
        return []


@pytest.mark.parametrize("chrom,start,end,expect", [
    ("chr1", 1, 100, (10, 110)),                 # 1-based inclusive in, 0-based half-open out, shifted
    ("chr1", 550, 700, None),                     # start in the hole: liftover drop
    ("chr2", 1, 100, None),                       # reversed block: liftover drop
    ("chr3", 1, 100, (0, 100)),                   # untouched
    ("chr3", 900_000, 1_100_000, "length"),       # 200 kb straddling a 5 Mb insertion: beyond 10% and 1 Mb
    ("chr3", 1, 20_000_000, "length"),            # 20 Mb gaining 5 Mb: 25% and 5 Mb, beyond both
    ("chr3", 1, 90_000_000, (0, 95_000_000)),     # 90 Mb gaining 5 Mb: 5.6%, within the fraction, kept
    ("chr3", 999_000, 1_000_500, (999_000 - 1, 6_000_500)),  # 1.5 kb gaining 5 Mb... beyond both: see next
])
def test_lift_segment_keeps_cleanly_lifted_segments_within_the_length_tolerance(chrom, start, end, expect):
    got = build.lift_segment(_Chain(), chrom, start, end)
    if (chrom, start) == ("chr3", 999_000):
        assert got == "length"                    # a tiny segment stretched to 5 Mb is the textbook case
    else:
        assert got == expect
