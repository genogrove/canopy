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
