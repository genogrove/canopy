# SPDX-License-Identifier: GPL-3.0-or-later
"""load_grove builds a .gg once, then deserializes it. Runs only where pygenogrove
is installed (CI skips); uses a file:// resource so it needs no network."""

from __future__ import annotations

import hashlib

import pytest

pg = pytest.importorskip("pygenogrove")

from genogrove_canopy import resources
from genogrove_canopy.resources import Resource

GFF3 = (
    "##gff-version 3\n"
    "chr1\tH\tgene\t1000\t2000\t.\t+\t.\tID=g1;gene_name=AAA\n"
    "chr1\tH\ttranscript\t1000\t2000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tH\texon\t1000\t1500\t.\t+\t.\tID=e1;Parent=t1\n"
).encode()


def _register(monkeypatch, tmp_path) -> str:
    src = tmp_path / "mini.gff3"
    src.write_bytes(GFF3)
    sha = hashlib.sha256(GFF3).hexdigest()
    monkeypatch.setattr(resources, "_CACHE", tmp_path / "cache")
    monkeypatch.setitem(resources.RESOURCES, "_mini", Resource("_mini", src.as_uri(), sha))
    return sha


def _all_gg(tmp_path, sha):
    return tmp_path / "cache" / "groves" / f"{sha}.{resources._GROVE_SCHEMA}" / "_all.gg"


def test_grove_index_builds_shards_and_whole(monkeypatch, tmp_path) -> None:
    sha = _register(monkeypatch, tmp_path)

    shards, all_path = resources.grove_index("_mini")  # first call builds the index
    assert set(shards) == {"chr1"}                      # the fixture is single-chromosome
    assert _all_gg(tmp_path, sha).exists()
    # the chr1 shard and the whole grove both hold the 3 features, edges intact
    chr1 = pg.Grove.deserialize(shards["chr1"])
    assert chr1.size() == 3
    gene = next(k for k in chr1.intersect(pg.GenomicCoordinate("*", 1500, 1500), "chr1")
                if k.data["type"] == "gene")
    assert [n.data["type"] for n in chr1.get_neighbors(gene)] == ["transcript"]
    assert pg.Grove.deserialize(all_path).size() == 3

    assert resources.grove_index("_mini")[1] == all_path  # second call: cache hit, no rebuild


def test_load_grove_returns_whole_and_self_heals(monkeypatch, tmp_path) -> None:
    sha = _register(monkeypatch, tmp_path)
    assert resources.load_grove("_mini").size() == 3      # whole-genome grove

    _all_gg(tmp_path, sha).write_bytes(b"not a grove")    # corrupt the cached index
    assert resources.load_grove("_mini").size() == 3      # self-heals: rebuild rather than crash
