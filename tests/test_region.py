# SPDX-License-Identifier: GPL-3.0-or-later
"""Region reads via build_grove on a bgzip+tabix GFF. Runs only where pygenogrove
(with region support) AND htslib's bgzip/tabix are available; skips otherwise."""

from __future__ import annotations

import shutil
import subprocess

import pytest

pg = pytest.importorskip("pygenogrove")
if not (shutil.which("bgzip") and shutil.which("tabix")):
    pytest.skip("bgzip/tabix (htslib) not on PATH", allow_module_level=True)

from genogrove_canopy.gff import build_grove

GFF = (
    "##gff-version 3\n"
    "chr1\tH\tgene\t1000\t2000\t.\t+\t.\tID=g1\n"
    "chr1\tH\tgene\t8000\t9000\t.\t-\t.\tID=g2\n"
)


def _bgzip_tabix(tmp_path):
    raw = tmp_path / "a.gff3"
    raw.write_text(GFF)
    gz = tmp_path / "a.gff3.gz"
    subprocess.run(
        f"(grep '^#' {raw}; grep -v '^#' {raw} | sort -k1,1 -k4,4n) | bgzip -c > {gz}",
        shell=True, check=True,
    )
    subprocess.run(["tabix", "-p", "gff", str(gz)], check=True)
    return str(gz)


def test_build_grove_reads_only_the_region(tmp_path):
    gz = _bgzip_tabix(tmp_path)

    # region "chr1:1500-1500" (1-based) overlaps g1 [1000,2000] but not g2 [8000,9000]
    g = build_grove(gz, region="chr1:1500-1500")
    genes = [k.data["id"] for k in g.intersect(pg.GenomicCoordinate("*", 1499, 1499), "chr1")
             if k.data["type"] == "gene"]
    assert genes == ["g1"]
    assert g.size() == 1                      # only the overlapping feature was loaded

    assert build_grove(gz).size() == 2        # region="" reads the whole file
