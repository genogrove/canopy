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
    assert "size 3" in result.stdout


def test_sandbox_still_blocks_network_with_pygenogrove_on_path(tmp_path):
    # The widened sys.path must not let the untrusted code import a blocked module.
    result = sandbox.run("import socket", extra_syspath=[_pygenogrove_site_dir()])
    assert result.returncode != 0
    assert "blocked" in result.stderr or "socket" in result.stderr
