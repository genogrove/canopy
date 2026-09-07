# SPDX-License-Identifier: GPL-3.0-or-later
"""SV layer: segment cutting, backbone anchoring, breakpoint edges. Runs only where
pygenogrove is installed (CI); skipped in the bare skeleton env."""

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


def _edge_neighbor(grove, key, sv_id):
    return next(iter(grove.get_neighbors_if(
        key, lambda m: m and m.get("rel") == "breakpoint_edge" and m.get("sv_id") == sv_id)))


def _anchor(grove, key):
    return next(iter(grove.get_neighbors_if(key, lambda m: m and m.get("rel") == "anchored_to")))


def test_translocation_anchors_both_genes_and_connects_segments():
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "A"})
    g.insert("chr2", pg.GenomicCoordinate("+", 5000, 6000), {"type": "gene", "id": "B"})

    records = [_sv("S1", "SV1", "chr1", 1500, "+", "chr2", 5500, "-", "TRA")]
    n, created = sv.attach_tracked(g, records)
    assert n == 1

    # "+" at 1500 keeps base 1500: the joined segment is chr1:[0, 1500]. "-" at 5500 continues
    # from base 5500: chr2:[5500, ...]. Each cut also leaves the far side as its own segment.
    segs = {(c, k.value.start, k.value.end) for c, k in created if k.data["type"] == "sv_segment"}
    assert segs == {("chr1", 0, 1500), ("chr1", 1501, 1502), ("chr2", 0, 5499), ("chr2", 5500, 5501)}
    seg1 = next(k for c, k in created if c == "chr1" and k.value.start == 0)
    seg2 = next(k for c, k in created if c == "chr2" and k.value.start == 5500)

    assert _anchor(g, seg1).data["id"] == "A"
    assert _anchor(g, seg2).data["id"] == "B"

    edge_target = _edge_neighbor(g, seg1, "SV1")
    assert edge_target.value.start == seg2.value.start
    edges = g.get_edges(seg1)
    match = next(e for e in edges if e.get("sv_id") == "SV1")
    assert match["svclass"] == "TRA" and match["orientation"] == "+/-"


def test_one_breakpoint_intergenic_creates_bin():
    g = pg.Grove(order=100)
    g.insert("chr3", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "C"})

    records = [_sv("S1", "SV1", "chr3", 1500, "+", "chr3", 50000, "-", "DEL")]
    n, created = sv.attach_tracked(g, records)
    assert n == 1

    bins = [k for _, k in created if k.data["type"] == "intergenic_region"]
    assert len(bins) == 1
    assert (bins[0].value.start, bins[0].value.end) == (0, 999_999)  # 0-based closed

    far_seg = next(k for _, k in created
                   if k.data["type"] == "sv_segment" and k.value.start == 50000)
    assert _anchor(g, far_seg).data["type"] == "intergenic_region"


def test_two_intergenic_breakpoints_share_one_bin():
    g = pg.Grove(order=100)
    records = [_sv("S1", "SV1", "chr4", 100, "+", "chr4", 200, "-", "DEL")]
    n, created = sv.attach_tracked(g, records)
    assert n == 1

    bins = [k for _, k in created if k.data["type"] == "intergenic_region"]
    assert len(bins) == 1  # both breakpoints land in bin 0 -> one bin, reused, not duplicated


def test_two_svs_sharing_one_gene_get_independent_edges():
    g = pg.Grove(order=100)
    g.insert("chr5", pg.GenomicCoordinate("+", 1000, 5000), {"type": "gene", "id": "D"})

    records = [
        _sv("S1", "SV1", "chr5", 1500, "+", "chr5", 20000, "-", "DEL"),
        _sv("S1", "SV2", "chr5", 3000, "+", "chr5", 25000, "-", "DUP"),
    ]
    n, created = sv.attach_tracked(g, records)
    assert n == 2

    # Cuts at 1501 and 3001: [0,1500] carries SV1's edge, [1501,3000] carries SV2's — both
    # overlap gene D and anchor to it independently, and [3001,19999] is the deleted middle.
    gene_segs = [k for _, k in created if k.data["type"] == "sv_segment"
                 and _anchor(g, k).data.get("id") == "D"]
    assert sorted((k.value.start, k.value.end) for k in gene_segs) == [(0, 1500), (1501, 3000), (3001, 19999)]
    sv_ids = {e["sv_id"] for k in gene_segs for e in g.get_edges(k) if e.get("rel") == "breakpoint_edge"}
    assert sv_ids == {"SV1", "SV2"}


def test_insertion_carries_payload_instead_of_second_segment():
    g = pg.Grove(order=100)
    records = [_sv("S1", "SV1", "chr6", 100, "+", "chr6", 101, "-", "INS",
                    length="50", insertion_class="MEI", sequence=None)]
    n, created = sv.attach_tracked(g, records)
    assert n == 1

    # 100(+) / 101(-) is one cut between bases 100 and 101: two segments, no orphan in between.
    segs = sorted((k.value.start, k.value.end) for _, k in created if k.data["type"] == "sv_segment")
    assert segs == [(0, 100), (101, 102)]
    seg1 = next(k for _, k in created if k.value.start == 0)
    match = next(e for e in g.get_edges(seg1) if e.get("sv_id") == "SV1")
    assert match["svclass"] == "INS"
    assert match["length"] == "50" and match["insertion_class"] == "MEI"
    assert match["sequence"] is None  # never invented when the caller couldn't resolve it


def test_detach_removes_only_what_attach_created():
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "A"})
    before = g.size()

    records = [_sv("S1", "SV1", "chr1", 1500, "+", "chr1", 50000, "-", "DEL")]
    n, created = sv.attach_tracked(g, records)
    assert g.size() > before

    sv.detach(g, created)
    assert g.size() == before
    hits = list(g.intersect(pg.GenomicCoordinate("*", 1000, 2000), "chr1"))
    assert len(hits) == 1 and hits[0].data["id"] == "A"  # the gene survives untouched


def test_parse_bedpe(tmp_path):
    p = tmp_path / "sample.bedpe"
    p.write_text(
        "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\t"
        "svclass\tsvmethod\n"
        "chr1\t1000\t1001\tchr1\t2000\t2001\tSV1\t5\t+\t-\tDEL\tm\n"
    )
    records = sv.parse_bedpe(p)
    assert records == [{
        "chrom1": "chr1", "start1": "1000", "end1": "1001",
        "chrom2": "chr1", "start2": "2000", "end2": "2001",
        "sv_id": "SV1", "pe_support": "5", "strand1": "+", "strand2": "-",
        "svclass": "DEL", "svmethod": "m",
    }]


def test_deletion_joins_flanks_and_orphans_the_deleted_middle():
    """The first '+' breakend on a chromosome used to fall through to the segment STARTING at the
    cut, so a lone DEL joined the deleted piece itself to the downstream flank. There is now a
    leading segment from 0, so upstream -> downstream is the edge and the middle is unreached."""
    g = pg.Grove(order=100)
    records = [_sv("S1", "SV1", "chr1", 1500, "+", "chr1", 20000, "-", "DEL")]
    _, created = sv.attach_tracked(g, records)
    by_start = {k.value.start: k for _, k in created if k.data["type"] == "sv_segment"}
    assert sorted((k.value.start, k.value.end) for k in by_start.values()) == [
        (0, 1500), (1501, 19999), (20000, 20001)]
    assert _edge_neighbor(g, by_start[0], "SV1").value.start == 20000
    assert not g.get_neighbors_if(by_start[1501], lambda m: m and m.get("rel") == "breakpoint_edge")
