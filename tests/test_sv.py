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

    seg1 = next(k for chrom, k in created if chrom == "chr1" and k.data["type"] == "sv_segment")
    seg2 = next(k for chrom, k in created if chrom == "chr2" and k.data["type"] == "sv_segment")

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
    assert (bins[0].value.start, bins[0].value.end) == (0, 1_000_000)

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

    gene_seg = next(k for _, k in created
                     if k.data["type"] == "sv_segment" and k.value.start == 1500)
    assert _anchor(g, gene_seg).data["id"] == "D"

    sv_ids = {e["sv_id"] for e in g.get_edges(gene_seg) if e.get("rel") == "breakpoint_edge"}
    assert sv_ids == {"SV1", "SV2"}  # both SVs reach the same gene-anchored segment, no conflict


def test_insertion_carries_payload_instead_of_second_segment():
    g = pg.Grove(order=100)
    records = [_sv("S1", "SV1", "chr6", 100, "+", "chr6", 101, "-", "INS",
                    length="50", insertion_class="MEI", sequence=None)]
    n, created = sv.attach_tracked(g, records)
    assert n == 1

    seg1 = next(k for _, k in created if k.value.start == 100)
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
