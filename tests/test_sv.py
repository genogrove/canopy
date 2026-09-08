# SPDX-License-Identifier: GPL-3.0-or-later
"""SV layer: breakend → gene-or-bin anchoring, one breakpoint edge per SV, exact detach. Runs
only where pygenogrove is installed (CI); skipped in the bare skeleton env."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from genogrove_canopy.layers import sv


def _sv(sample, sv_id, chrom1, start1, strand1, chrom2, start2, strand2, svclass, **extra):
    rec = {
        "chrom1": chrom1, "start1": str(start1), "end1": str(start1 + 1),
        "chrom2": chrom2, "start2": str(start2), "end2": str(start2 + 1),
        "sv_id": sv_id, "pe_support": "10", "strand1": strand1, "strand2": strand2,
        "svclass": svclass, "svmethod": "m", "sample": sample,
    }
    rec.update(extra)
    return rec


def _bp_edges(grove, key):
    """(target, payload) for every breakpoint edge leaving ``key``."""
    return [(t, m) for t, m in grove.get_edge_list(key) if m and m.get("rel") == "breakpoint_edge"]


def _gene(grove, chrom, pos):
    return next(k for k in grove.intersect(pg.GenomicCoordinate("*", pos, pos), chrom)
                if k.data.get("type") == "gene")


def test_translocation_joins_the_two_genes_with_positions_on_the_edge():
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "A"})
    g.insert("chr2", pg.GenomicCoordinate("+", 5000, 6000), {"type": "gene", "id": "B"})

    n, created = sv.attach_tracked(g, [_sv("S1", "SV1", "chr1", 1500, "+", "chr2", 5500, "-", "TRA")])
    assert n == 1
    a, b = _gene(g, "chr1", 1500), _gene(g, "chr2", 5500)
    (t, m), = _bp_edges(g, a)
    assert t.data["id"] == "B" and _bp_edges(g, b)[0][0].data["id"] == "A"   # both directions
    assert m["svclass"] == "TRA" and m["sample"] == "S1"
    assert (m["chrom1"], m["pos1"], m["strand1"]) == ("chr1", 1500, "+")
    assert (m["chrom2"], m["pos2"], m["strand2"]) == ("chr2", 5500, "-")
    assert created == [("edge", a, b)]                    # nothing inserted: genes are the anchors


def test_breakend_inside_overlapping_genes_anchors_to_each():
    g = pg.Grove(order=100)
    g.insert("chr7", pg.GenomicCoordinate("+", 1000, 3000), {"type": "gene", "id": "EGFR"})
    g.insert("chr7", pg.GenomicCoordinate("-", 2500, 3500), {"type": "gene", "id": "EGFR-AS1"})
    g.insert("chr9", pg.GenomicCoordinate("+", 100, 200), {"type": "gene", "id": "C"})

    _, created = sv.attach_tracked(g, [_sv("S1", "SV1", "chr7", 2700, "+", "chr9", 150, "-", "TRA")])
    c = _gene(g, "chr9", 150)
    assert sorted(t.data["id"] for t, _ in _bp_edges(g, c)) == ["EGFR", "EGFR-AS1"]
    assert len(created) == 2 and all(e[0] == "edge" for e in created)


def test_intergenic_breakend_creates_a_closed_1mb_bin():
    g = pg.Grove(order=100)
    g.insert("chr3", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "C"})

    _, created = sv.attach_tracked(g, [_sv("S1", "SV1", "chr3", 1500, "+", "chr3", 2_500_000, "-", "DEL")])
    bins = [k for e, _, k in [c for c in created if c[0] == "key"] if k.data["type"] == "intergenic_region"]
    assert len(bins) == 1
    assert (bins[0].value.start, bins[0].value.end) == (2_000_000, 2_999_999)  # 0-based closed
    (t, m), = _bp_edges(g, _gene(g, "chr3", 1500))
    assert t.data["type"] == "intergenic_region" and m["pos2"] == 2_500_000


def test_two_intergenic_breakends_share_one_bin():
    g = pg.Grove(order=100)
    _, created = sv.attach_tracked(g, [_sv("S1", "SV1", "chr4", 100, "+", "chr4", 200, "-", "DEL")])
    bins = [c for c in created if c[0] == "key"]
    assert len(bins) == 1                                 # same bin at both ends: reused, not duplicated
    (_, _, b), = bins
    assert [(t.value.start, m["sv_id"]) for t, m in _bp_edges(g, b)] == [(0, "SV1")]  # one self-edge


def test_two_svs_on_one_gene_are_two_edges():
    g = pg.Grove(order=100)
    g.insert("chr5", pg.GenomicCoordinate("+", 1000, 5000), {"type": "gene", "id": "D"})
    g.insert("chr5", pg.GenomicCoordinate("+", 20_000, 21_000), {"type": "gene", "id": "E"})

    n, _ = sv.attach_tracked(g, [
        _sv("S1", "SV1", "chr5", 1500, "+", "chr5", 20_500, "-", "DEL"),
        _sv("S1", "SV2", "chr5", 3000, "+", "chr5", 20_700, "-", "DUP"),
    ])
    assert n == 2
    edges = _bp_edges(g, _gene(g, "chr5", 1500))
    assert sorted((m["sv_id"], m["svclass"], m["pos1"]) for _, m in edges) == [
        ("SV1", "DEL", 1500), ("SV2", "DUP", 3000)]


def test_insertion_carries_its_payload_on_the_edge():
    g = pg.Grove(order=100)
    g.insert("chr6", pg.GenomicCoordinate("+", 50, 150), {"type": "gene", "id": "F"})
    sv.attach_tracked(g, [_sv("S1", "SV1", "chr6", 100, "+", "chr6", 101, "-", "INS",
                               length="50", insertion_class="MEI", sequence=None)])
    (t, m), = _bp_edges(g, _gene(g, "chr6", 100))          # both breakends in F: ONE self-edge
    assert t.data["id"] == "F"
    assert m["svclass"] == "INS" and m["length"] == "50" and m["insertion_class"] == "MEI"
    assert m["sequence"] is None                          # never invented


def test_detach_removes_only_what_attach_created():
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "A", "name": "A"})
    g.insert("chr1", pg.GenomicCoordinate("+", 3000, 4000), {"type": "gene", "id": "B"})
    tx = g.add_external_key(pg.GenomicCoordinate("+", 1000, 2000), {"type": "transcript"})
    g.add_edge(_gene(g, "chr1", 1500), tx, {"rel": "contains"})          # pre-existing backbone edge
    before = (g.size(), g.vertex_count(), g.edge_count())

    _, created = sv.attach_tracked(g, [
        _sv("S1", "SV1", "chr1", 1500, "+", "chr1", 5_500_000, "-", "DEL"),   # gene -> new bin
        _sv("S1", "SV2", "chr1", 1600, "+", "chr1", 3500, "-", "DUP"),        # gene -> gene
    ])
    assert g.edge_count() == before[2] + 4 and g.size() == before[0] + 1

    sv.detach(g, created)
    assert (g.size(), g.vertex_count(), g.edge_count()) == before
    a = _gene(g, "chr1", 1500)
    assert a.data["name"] == "A" and [m for _, m in g.get_edge_list(a)] == [{"rel": "contains"}]


def test_parse_bedpe(tmp_path):
    p = tmp_path / "sample.bedpe"
    p.write_text(
        "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\t"
        "svclass\tsvmethod\n"
        "chr1\t1000\t1001\tchr1\t2000\t2001\tSV1\t5\t+\t-\tDEL\tm\n"
    )
    assert sv.parse_bedpe(p) == [{
        "chrom1": "chr1", "start1": "1000", "end1": "1001",
        "chrom2": "chr1", "start2": "2000", "end2": "2001",
        "sv_id": "SV1", "pe_support": "5", "strand1": "+", "strand2": "-",
        "svclass": "DEL", "svmethod": "m",
    }]
