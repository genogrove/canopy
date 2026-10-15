# SPDX-License-Identifier: GPL-3.0-or-later
"""Expression layer: one graph-only node per gene per cohort hung off the gene, values keyed by
the SV layer's sample id, exact detach, cohort extraction. Runs only where pygenogrove is
installed."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from genogrove_canopy.layers import expression


def _backbone():
    g = pg.Grove(order=100)
    g.insert("chr5", pg.GenomicCoordinate("-", 1_253_147, 1_295_068), {"type": "gene", "id": "ENSG00000164362.21", "name": "TERT"})
    g.insert("chr8", pg.GenomicCoordinate("+", 127_735_433, 127_742_951), {"type": "gene", "id": "ENSG00000136997.20", "name": "MYC"})
    return g


def _table(tmp_path, name="LIRI-JP.tsv"):
    p = tmp_path / name
    # TSS is 1-based as in gene_tss.tsv.gz: TERT (minus strand) ends at 0-based 1,295,068, so its
    # 1-based TSS is 1,295,069 — one base outside the gene without the -1 in attach_tracked
    p.write_text("gene\tchrom\ttss\tA1\tA2\tA3\n"
                 "ENSG00000164362\tchr5\t1295069\t6.1\t-9.9\tNA\n"     # TERT: minus strand, TSS = gene end
                 "ENSG00000136997\tchr8\t127735434\t8.0\t7.5\t7.9\n"   # MYC: plus strand, TSS = gene start
                 "ENSG00000000000\tchr1\t100\t1\t2\t3\n")              # not in this backbone
    return p


def _expr(g, gene):
    return [(t, m) for t, m in g.get_edge_list(gene) if m and m.get("rel") == "expressed"]


def test_values_hang_off_the_gene_keyed_by_sample_and_never_come_back_from_intersect(tmp_path):
    g = _backbone()
    n, created, missed = expression.attach_tracked(g, _table(tmp_path), "LIRI-JP")
    assert (n, missed) == (2, 1) and [c[0] for c in created] == ["node", "node"]
    assert g.external_vertex_count() == 2 and g.edge_count() == 4
    tert, = [k for k in g.intersect(pg.GenomicCoordinate("*", 1_260_000, 1_260_000), "chr5") if k.data.get("type") == "gene"]
    (node, m), = _expr(g, tert)
    assert m["cohort"] == "LIRI-JP" and node.data["gene"] == "ENSG00000164362"
    assert node.data["values"] == {"A1": 6.1, "A2": -9.9} and node.data["unit"] == "log2(fpkm-uq+0.001)"   # NA dropped
    back, = [(t, mm) for t, mm in g.get_edge_list(node) if mm and mm["rel"] == "expression_of"]
    assert back[0].data["name"] == "TERT"
    assert all(k.data["type"] == "gene" for k in g.intersect(pg.GenomicCoordinate("*", 1_253_147, 1_295_068), "chr5"))
    # a second cohort is its own node on the same gene
    expression.attach_tracked(g, _table(tmp_path, "LINC-JP.tsv"), "LINC-JP")
    assert sorted(m["cohort"] for _, m in _expr(g, tert)) == ["LINC-JP", "LIRI-JP"]


def test_detach_removes_this_cohorts_edges_only(tmp_path):
    g = _backbone()
    _, c1, _ = expression.attach_tracked(g, _table(tmp_path), "LIRI-JP")
    _, c2, _ = expression.attach_tracked(g, _table(tmp_path, "LINC-JP.tsv"), "LINC-JP")
    expression.detach(g, c1)
    tert, = [k for k in g.intersect(pg.GenomicCoordinate("*", 1_260_000, 1_260_000), "chr5") if k.data.get("type") == "gene"]
    assert [m["cohort"] for _, m in _expr(g, tert)] == ["LINC-JP"]
    expression.detach(g, c2)
    assert g.edge_count() == 0


def test_cohort_file_extracts_one_project_from_the_pinned_tarball(monkeypatch, tmp_path):
    import tarfile

    from genogrove_canopy import resources

    table = _table(tmp_path)
    tgz = tmp_path / "expr.tgz"
    with tarfile.open(tgz, "w:gz") as t:
        t.add(table, arcname="LIRI-JP.tsv")
    monkeypatch.setattr(resources, "resolve", lambda name: {"pcawg.expression": tgz}[name])
    monkeypatch.setattr(resources, "expression_cohorts", lambda: [{"project_code": "LIRI-JP", "n_samples": 3, "n_genes": 2}])
    monkeypatch.setattr(expression, "EXPRESSION_DIR", tmp_path / "expression_cohorts")
    out = expression.cohort_file("LIRI-JP")
    assert out.parent.parent == (tmp_path / "expression_cohorts").resolve() and out.name == "LIRI-JP.tsv"
    assert out.read_text().splitlines()[0] == "gene\tchrom\ttss\tA1\tA2\tA3"
    with pytest.raises(KeyError, match="not a PCAWG expression cohort"):
        expression.cohort_file("NOPE-XX")


def test_a_re_pin_never_serves_a_table_extracted_under_the_old_pin(monkeypatch, tmp_path):
    import dataclasses
    import tarfile

    from genogrove_canopy import resources

    tgz = tmp_path / "expr.tgz"
    with tarfile.open(tgz, "w:gz") as t:
        t.add(_table(tmp_path), arcname="LIRI-JP.tsv")
    monkeypatch.setattr(resources, "resolve", lambda name: {"pcawg.expression": tgz}[name])
    monkeypatch.setattr(resources, "expression_cohorts", lambda: [{"project_code": "LIRI-JP", "n_samples": 3, "n_genes": 2}])
    monkeypatch.setattr(expression, "EXPRESSION_DIR", tmp_path / "expression_cohorts")
    first = expression.cohort_file("LIRI-JP")
    first.write_text("stale")                      # what a user extracted under the old pin
    assert expression.cohort_file("LIRI-JP").read_text() == "stale"   # same pin: cached
    old = resources.RESOURCES["pcawg.expression"]
    monkeypatch.setitem(resources.RESOURCES, "pcawg.expression", dataclasses.replace(old, sha256="f" * 64))
    second = expression.cohort_file("LIRI-JP")
    assert second != first and second.read_text().startswith("gene\tchrom\ttss")   # new pin: re-extracted
