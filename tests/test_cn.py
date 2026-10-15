# SPDX-License-Identifier: GPL-3.0-or-later
"""Copy-number layer: one indexed segment per sample, reached by intersect and filtered on
sample, null allele values kept null, exact detach, cohort extraction. Runs only where
pygenogrove is installed."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from genogrove_canopy.layers import cn


def _seg(chrom, start, end, sample, total, major="", minor="", star="", cohort="BRCA-US"):
    return {"chrom": chrom, "start": str(start), "end": str(end), "total_cn": str(total),
            "major_cn": str(major), "minor_cn": str(minor), "star": str(star), "sample": sample, "cohort": cohort}


def _backbone():
    g = pg.Grove(order=100)
    g.insert("chr8", pg.GenomicCoordinate("+", 127_735_433, 127_742_951), {"type": "gene", "id": "MYC", "name": "MYC"})
    return g


def _under(g, start, end):
    return [k for k in g.intersect(pg.GenomicCoordinate("*", start, end), "chr8") if k.data.get("type") == "cn_segment"]


def test_segments_are_intervals_filtered_by_sample_with_breakpoint_genes_showing_both():
    g = _backbone()
    n, created = cn.attach_tracked(g, [
        _seg("chr8", 127_000_000, 128_000_000, "S1", 6, 5, 1, 3),      # S1: MYC amplified, one segment
        _seg("chr8", 127_000_000, 127_740_000, "S2", 2, 1, 1, 2),      # S2: breakpoint inside MYC
        _seg("chr8", 127_740_000, 128_000_000, "S2", 3, 3, 0, 2),      #     -> two segments under MYC
        _seg("chr8", 100, 200, "S3", 1),                                 # elsewhere, no allele values
        _seg("chr8", 300, 400, "S3", "NA"),                              # the consensus made no call
    ])
    assert n == 5 and [c[0] for c in created] == ["key"] * 5 and g.external_vertex_count() == 0
    s3_uncalled, = _under(g, 350, 350)
    assert s3_uncalled.data["total_cn"] is None
    assert g.edge_count() == 0                                                  # overlap is spatial only
    under = _under(g, 127_735_433, 127_742_951)
    by_sample = {}
    for k in under:
        by_sample.setdefault(k.data["sample"], []).append(k.data["total_cn"])
    assert by_sample == {"S1": [6], "S2": [2, 3]}
    s1, = [k for k in under if k.data["sample"] == "S1"]
    assert (s1.value.start, s1.value.end) == (127_000_000, 127_999_999)        # BED half-open -> closed
    assert (s1.data["major_cn"], s1.data["minor_cn"], s1.data["star"], s1.data["source"]) == (5, 1, 3, "PCAWG-consensus-CN")
    s3, = _under(g, 150, 150)
    assert (s3.data["major_cn"], s3.data["minor_cn"], s3.data["star"]) == (None, None, None)   # null, never 0


def test_detach_removes_exactly_this_attachments_segments():
    g = _backbone()
    _, c1 = cn.attach_tracked(g, [_seg("chr8", 127_000_000, 128_000_000, "S1", 6)])
    _, c2 = cn.attach_tracked(g, [_seg("chr8", 127_000_000, 128_000_000, "S2", 2)])
    cn.detach(g, c1)
    assert [k.data["sample"] for k in _under(g, 127_740_000, 127_740_000)] == ["S2"]
    cn.detach(g, c2)
    assert _under(g, 127_740_000, 127_740_000) == [] and g.vertex_count() == 1


def test_cohort_file_extracts_one_project_from_the_pinned_tarball(monkeypatch, tmp_path):
    import tarfile

    from genogrove_canopy import resources

    (tmp_path / "BRCA-US.tsv").write_text("\t".join(cn._FIELDS) + "\nchr8\t0\t1000\t2\t1\t1\t3\tS1\tBRCA-US\n")
    tgz = tmp_path / "cn.tgz"
    with tarfile.open(tgz, "w:gz") as t:
        t.add(tmp_path / "BRCA-US.tsv", arcname="BRCA-US.tsv")
    monkeypatch.setattr(resources, "resolve", lambda name: {"pcawg.cn": tgz}[name])
    monkeypatch.setattr(resources, "cn_cohorts", lambda: [{"project_code": "BRCA-US", "n_samples": 1, "n_segments": 1,
                                                            "n_dropped_liftover": 0, "n_dropped_length": 0}])
    monkeypatch.setattr(cn, "CN_DIR", tmp_path / "cn_cohorts")
    out = cn.cohort_file("BRCA-US")
    assert out.parent.parent == (tmp_path / "cn_cohorts").resolve() and out.name == "BRCA-US.tsv"
    assert out.read_text().splitlines()[1].endswith("\tS1\tBRCA-US")
    with pytest.raises(KeyError, match="not a PCAWG copy-number cohort"):
        cn.cohort_file("NOPE-XX")
