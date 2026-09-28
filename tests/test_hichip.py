# SPDX-License-Identifier: GPL-3.0-or-later
"""HiChIP layer: windows as indexed keys inserted on demand, no overlap edges, one contact
edge per loop per sample, exact detach. Runs only where pygenogrove is installed."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from genogrove_canopy.layers import hichip


def _loop(sample, chrom, s1, s2, cc=10, q=0.01, **extra):
    rec = {"chrom1": chrom, "start1": str(s1), "end1": str(s1 + 10_000),
           "chrom2": chrom, "start2": str(s2), "end2": str(s2 + 10_000),
           "cc": str(cc), "q": str(q), "peak1": "1", "peak2": "0", "sample": sample, "cohort": "BRCA"}
    rec.update(extra)
    return rec


def _span(g, chrom, start, end, type_):
    return [k for k in g.intersect(pg.GenomicCoordinate("*", start, end), chrom) if k.data.get("type") == type_]


def _contacts(grove, window):
    return [(t, m) for t, m in grove.get_edge_list(window) if m and m.get("rel") == "contact_edge"]


def _backbone():
    g = pg.Grove(order=100)
    g.insert("chr8", pg.GenomicCoordinate("+", 127_735_000, 127_742_000), {"type": "gene", "id": "MYC"})
    g.insert("chr8", pg.GenomicCoordinate(".", 127_731_000, 127_731_300), {"type": "regulatory_region", "id": "E1"})
    g.insert("chr8", pg.GenomicCoordinate("+", 127_900_000, 128_200_000), {"type": "gene", "id": "PVT1"})
    return g


def test_loop_is_reachable_from_gene_by_intersect_through_window_to_partner():
    g = _backbone()
    before = g.vertex_count()
    n, created = hichip.attach_tracked(g, [_loop("S1", "chr8", 127_730_000, 128_100_000)])
    assert n == 1
    assert g.vertex_count() == before + 2 and g.external_vertex_count() == 0   # windows are indexed
    assert g.edge_count() == 2                                                 # the loop only, no overlap edges

    myc, = _span(g, "chr8", 127_735_000, 127_742_000, "gene")
    w, = _span(g, "chr8", myc.value.start, myc.value.end, "hichip_anchor")     # gene -> window: spatial
    assert (w.value.start, w.value.end, w.data["chrom"]) == (127_730_000, 127_739_999, "chr8")  # 0-based closed
    (partner, m), = _contacts(g, w)
    assert m["sample"] == "S1" and m["cc"] == 10 and m["q"] == 0.01 and (m["peak1"], m["peak2"]) == (1, 0)
    assert _contacts(g, partner)[0][0] is w                                     # both directions
    other_end = g.intersect(  # "*": a "." query misses stranded genes
        pg.GenomicCoordinate("*", partner.value.start, partner.value.end), partner.data["chrom"])
    assert sorted(k.data["type"] for k in other_end) == ["gene", "hichip_anchor"]  # itself + PVT1
    assert [k.data["id"] for k in other_end if k.data["type"] == "gene"] == ["PVT1"]  # window inside a big gene
    assert _span(g, "chr8", 127_731_100, 127_731_100, "regulatory_region") and \
        _span(g, "chr8", 127_731_100, 127_731_100, "hichip_anchor") == [w]     # cCRE and window co-located, unlinked
    assert [e[0] for e in created] == ["key", "key", "edge"]


def test_window_is_shared_across_samples_and_across_caches():
    g = _backbone()
    windows = {}
    hichip.attach_tracked(g, [_loop("S1", "chr8", 127_730_000, 128_100_000)], windows)
    hichip.attach_tracked(g, [_loop("S2", "chr8", 127_730_000, 128_100_000, cc=3)], windows)
    _, created = hichip.attach_tracked(g, [_loop("S3", "chr8", 127_730_000, 127_900_000)])  # fresh cache
    assert len(windows) == 2 and g.vertex_count() == 3 + 3                     # found by intersect, not re-inserted
    assert [e[0] for e in created] == ["key", "edge"]                          # only the new window inserted
    w, = _span(g, "chr8", 127_735_000, 127_735_000, "hichip_anchor")
    assert sorted(m["sample"] for _, m in _contacts(g, w)) == ["S1", "S2", "S3"]


def test_cohort_file_extracts_one_project_from_the_pinned_tarball(monkeypatch, tmp_path):
    import tarfile

    from genogrove_canopy import resources

    (tmp_path / "TCGA-KIRC.tsv").write_text("\t".join(hichip._FIELDS) + "\nchr1\t0\t10000\tchr1\t50000\t60000\t5\t0.01\t1\t0\tS1\tTCGA-KIRC\n")
    tgz = tmp_path / "loops.tgz"
    with tarfile.open(tgz, "w:gz") as t:
        t.add(tmp_path / "TCGA-KIRC.tsv", arcname="TCGA-KIRC.tsv")
    monkeypatch.setattr(resources, "resolve", lambda name: {"tcga.hichip.loops": tgz}[name])
    monkeypatch.setattr(hichip, "HICHIP_DIR", tmp_path / "hichip_cohorts")
    out = hichip.cohort_file("TCGA-KIRC")
    assert out == (tmp_path / "hichip_cohorts" / "TCGA-KIRC.tsv").resolve()
    assert out.read_text().splitlines()[1].endswith("\tS1\tTCGA-KIRC")
    with pytest.raises(KeyError, match="not a TCGA HiChIP cohort"):
        hichip.cohort_file("TCGA-NOPE")


def test_detach_removes_only_this_attachment_and_its_orphaned_windows():
    g = _backbone()
    windows = {}
    _, c1 = hichip.attach_tracked(g, [_loop("S1", "chr8", 127_730_000, 128_100_000)], windows)
    _, c2 = hichip.attach_tracked(g, [_loop("S2", "chr8", 127_730_000, 127_900_000)], windows)
    hichip.detach(g, c1, windows)
    w, = _span(g, "chr8", 127_735_000, 127_735_000, "hichip_anchor")           # shared window kept
    assert [m["sample"] for _, m in _contacts(g, w)] == ["S2"]
    assert _span(g, "chr8", 128_100_000, 128_100_000, "hichip_anchor") == []   # S1-only window removed
    assert set(windows) == {("chr8", 127_730_000), ("chr8", 127_900_000)}
    hichip.detach(g, c2, windows)
    assert g.edge_count() == 0 and windows == {} and g.vertex_count() == 3
