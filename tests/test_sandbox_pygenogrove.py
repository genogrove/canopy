# SPDX-License-Identifier: GPL-3.0-or-later
"""The sandbox must actually run pygenogrove: import it (despite -S) and
deserialize a .gg. Validates the extra_syspath fix + the host/agent .gg path.
Runs only where pygenogrove is installed (CI skips)."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from genogrove_canopy import sandbox
from genogrove_canopy.cli import _pygenogrove_site_dir
from genogrove_canopy.gff import load_gff

GFF3 = (
    "##gff-version 3\n"
    "chr1\tH\tgene\t1000\t2000\t.\t+\t.\tID=g1;gene_name=AAA\n"
    "chr1\tH\ttranscript\t1000\t2000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tH\texon\t1000\t1500\t.\t+\t.\tID=e1;Parent=t1\n"
).encode()


def test_sandbox_deserializes_and_queries_a_grove(tmp_path):
    src = tmp_path / "mini.gff3"
    src.write_bytes(GFF3)
    gg = tmp_path / "mini.gg"
    load_gff(src).serialize(str(gg))

    # This is the shape of what the CLI runs: a host-injected path var + agent code.
    code = (
        f"GG = {str(gg)!r}\n"
        "import pygenogrove as pg\n"
        "g = pg.Grove.deserialize(GG)\n"
        "q = pg.GenomicCoordinate('*', 1200, 1200)\n"
        "print('genes', sum(1 for k in g.intersect(q, 'chr1') if k.data['type'] == 'gene'))\n"
        "print('size', g.size())\n"
    )
    result = sandbox.run(
        code, data_paths={"g": str(gg)}, extra_syspath=[_pygenogrove_site_dir()]
    )

    assert result.returncode == 0, result.stderr
    assert "genes 1" in result.stdout
    assert "size 1" in result.stdout  # only the gene is indexed; transcript + exon external


def test_sandbox_still_blocks_network_with_pygenogrove_on_path(tmp_path):
    # The widened sys.path must not let the untrusted code import a blocked module.
    result = sandbox.run("import socket", extra_syspath=[_pygenogrove_site_dir()])
    assert result.returncode != 0
    assert "blocked" in result.stderr or "socket" in result.stderr


def test_enhancer_preamble_attaches_links_inside_the_sandbox(tmp_path, monkeypatch):
    """The cohort preamble opens a links table from `enhancers.LINKS_DIR` inside the sandbox,
    so `_grove_context` must grant that directory — with only the grove granted, every enhancer
    question died on `PermissionError` and no test executed the preamble to notice."""
    from genogrove_canopy import resources
    from genogrove_canopy.cli import _grove_context
    from genogrove_canopy.layers import enhancers

    gg = tmp_path / "mini.gg"
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 999, 1999), {"type": "gene", "id": "ENSG1.2"})
    g.serialize(str(gg))
    links = tmp_path / "links" / "c.links.tsv"
    links.parent.mkdir()
    links.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t1\t0.5\t0.9\n")
    monkeypatch.setattr(enhancers, "LINKS_DIR", links.parent)
    monkeypatch.setattr(resources, "ensure_all_grove", lambda name: gg)

    _, _, data_paths = _grove_context()
    code = enhancers.preamble(str(gg), {"COH": str(links)}) + (
        "enh = [k for k in GROVE.intersect(pg.GenomicCoordinate('*', 150, 150), 'chr1')]\n"
        "print('enh', len(enh), enh[0].data['source'])\n"
    )
    result = sandbox.run(code, data_paths=data_paths, extra_syspath=[_pygenogrove_site_dir()])

    assert result.returncode == 0, result.stderr
    assert "enh 1 ENCODE-rE2G" in result.stdout


def test_warm_worker_reuses_the_attached_grove_and_a_query_cannot_mutate_it(tmp_path):
    """The path `serve` and `-i` always take: the cohort preamble runs twice in one `Worker`.
    The second run must hit `_CANOPY_STATE` (no re-attach) and see exactly what the first saw —
    even though the first query tried to insert into `GROVE`, which the read-only view refuses."""
    from genogrove_canopy.layers import enhancers

    gg = tmp_path / "mini.gg"
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 999, 1999), {"type": "gene", "id": "ENSG1.2"})
    g.serialize(str(gg))
    links = tmp_path / "c.links.tsv"
    links.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t1\t0.5\t0.9\n")
    pre = "import json\n" + enhancers.preamble(str(gg), {"C": str(links)})
    count = ("enh = [k for k in GROVE.intersect(pg.GenomicCoordinate('*', 0, 5000), 'chr1')"
             " if k.data.get('type') == 'enhancer']\nprint('enh', len(enh))\n")

    # The view hides the grove, but a bound method's `__self__` still names it — good enough to
    # prove both runs used the same object, i.e. the second was a memo hit.
    ident = "print('id', id(GROVE.intersect.__self__))\n"
    links2 = tmp_path / "d.links.tsv"       # a second cohort: the same element plus a new one
    links2.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t2\t0.4\t0.7\n"
                      "chr1\t300\t400\tintergenic\tENSG1\tchr1\t1000\t1\t0.1\t0.2\n")
    pre2 = "import json\n" + enhancers.preamble(str(gg), {"D": str(links2)})
    both = ("gene = next(iter(GROVE.intersect(pg.GenomicCoordinate('*', 1500, 1500), 'chr1')))\n"
            "print('cohorts', sorted(set(c for m in GROVE.get_edges(gene)"
            " if m and m.get('rel') == 'regulated_by' for c in m['byCohort'])), COHORTS)\n")
    w = sandbox.Worker(data_paths=[str(tmp_path)], extra_syspath=[_pygenogrove_site_dir()])
    try:
        first = w.submit(pre + count + ident
                         + "try:\n    GROVE.insert('chr1', pg.GenomicCoordinate('.', 5, 6), {})\n"
                         "except AttributeError as e:\n    print('refused', e)\n")
        second = w.submit(pre + count + ident)
        third = w.submit(pre2 + count + ident + both)
    finally:
        w.close()

    assert first.returncode == 0, first.stderr
    assert "enh 1" in first.stdout and "refused" in first.stdout
    assert second.returncode == 0, second.stderr
    assert "enh 1" in second.stdout                       # the refused insert left no trace
    assert third.returncode == 0, third.stderr
    # Cohort D attached onto the SAME grove: its new element appears, the shared element merged
    # (one node, both cohorts on the gene's edges), and COHORTS names only this question's.
    assert "enh 2" in third.stdout
    assert "cohorts ['C', 'D'] ['D']" in third.stdout
    ids = [ln for ln in (first.stdout + second.stdout + third.stdout).splitlines() if ln.startswith("id ")]
    assert len(ids) == 3 and len(set(ids)) == 1            # same grove object throughout
