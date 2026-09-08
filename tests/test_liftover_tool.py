# SPDX-License-Identifier: GPL-3.0-or-later
"""The one-off liftover tool's per-record logic, with a stub lifter (no pyliftover, no chain)."""
import gzip
import importlib.util
import io
import tarfile
from pathlib import Path

import pytest

_spec = importlib.util.spec_from_file_location(
    "liftover_pcawg_sv", Path(__file__).parents[1] / "tools" / "liftover_pcawg_sv.py")
lift = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lift)


class _Lifter:
    """chrom, pos -> [(chrom, pos + 10, strand, score)]; pos 999 maps nowhere, pos 500 reverses."""
    def convert_coordinate(self, chrom, pos):
        if pos == 999:
            return []
        return [(chrom, pos + 10, "-" if pos == 500 else "+", 1.0)]


def test_lift_one_reports_position_and_reversal():
    lo = _Lifter()
    assert lift._lift_one(lo, "1", 100) == (110, False)     # "1" -> "chr1", shifted, same strand
    assert lift._lift_one(lo, "chr1", 500) == (510, True)   # minus-strand chain block
    assert lift._lift_one(lo, "1", 999) is None             # no mapping


def test_lift_tarball_flips_the_strand_of_a_reversed_breakend(tmp_path):
    header = "\t".join(lift._FIELDS) + "\n"
    rows = ("1\t100\t101\t2\t500\t501\tSV1\t5\t+\t-\tTRA\tm\n"    # breakend 2 reverses
            "1\t100\t101\t1\t999\t1000\tSV2\t5\t+\t-\tDEL\tm\n")  # breakend 2 unmappable
    src = tmp_path / "in.tgz"
    with tarfile.open(src, "w:gz") as t:
        data = gzip.compress((header + rows).encode())
        info = tarfile.TarInfo("icgc/open/s.bedpe.gz"); info.size = len(data)
        t.addfile(info, io.BytesIO(data))

    out = tmp_path / "out.tgz"
    assert lift.lift_tarball(src, out, _Lifter()) == (1, 1)
    with tarfile.open(out) as t:
        text = gzip.decompress(t.extractfile("icgc/open/s.bedpe.gz").read()).decode()
    assert text.splitlines()[1] == "chr1\t110\t111\tchr2\t510\t511\tSV1\t5\t+\t+\tTRA\tm"


def test_chain_file_is_verified_even_when_cached(tmp_path, monkeypatch):
    import hashlib
    import pytest

    good = b"chain-bytes"
    monkeypatch.setattr(lift, "CHAIN_SHA256", hashlib.sha256(good).hexdigest())
    monkeypatch.setattr(lift.urllib.request, "urlretrieve", lambda url, dst: Path(dst).write_bytes(good))

    dest = tmp_path / "hg19ToHg38.over.chain.gz"
    dest.write_bytes(good[:4])                       # a truncated cached copy from an interrupted run
    with pytest.raises(RuntimeError, match="sha256 mismatch"):
        lift._chain_file(tmp_path)
    assert not dest.exists()                         # the bad copy is gone, not trusted next time

    assert lift._chain_file(tmp_path) == dest        # fresh download: verified, published
    assert dest.read_bytes() == good and not list(tmp_path.glob("*.part"))
    assert lift._chain_file(tmp_path) == dest        # cached good copy: verified again, kept


@pytest.mark.parametrize("source_class,strands,hits,expected_class,expected_junction", [
    # One breakend reverses: the source inversion call no longer has inversion geometry.
    ("t2tINV", ("-", "-"), ((110, "+"), (510, "-")), "INV", "DUP-like"),
    ("h2hINV", ("+", "+"), ((110, "-"), (510, "+")), "INV", "DUP-like"),
    # A reversed region swaps coordinate order as well as both strands.
    ("DEL", ("+", "-"), ((900, "-"), (500, "-")), "DEL", "DEL-like"),
    ("h2hINV", ("+", "+"), ((900, "-"), (500, "-")), "INV", "t2tINV"),
])
def test_lifted_call_keeps_provenance_separate_from_junction_geometry(
    tmp_path, source_class, strands, hits, expected_class, expected_junction,
):
    pg = pytest.importorskip("pygenogrove")
    from genogrove_canopy.layers import sv

    class Chain:
        def convert_coordinate(self, chrom, pos):
            mapped, strand = hits[0 if pos == 100 else 1]
            return [(chrom, mapped, strand, 1.0)]

    row = ["1", "100", "101", "1", "500", "501", "SV1", "5",
           *strands, source_class, "m"]
    src, out = tmp_path / "source.tgz", tmp_path / "lifted.tgz"
    with tarfile.open(src, "w:gz") as t:
        data = gzip.compress(("\t".join(lift._FIELDS) + "\n" + "\t".join(row) + "\n").encode())
        info = tarfile.TarInfo("sample.bedpe.gz")
        info.size = len(data)
        t.addfile(info, io.BytesIO(data))
    assert lift.lift_tarball(src, out, Chain()) == (1, 0)
    bedpe = tmp_path / "sample.bedpe.gz"
    with tarfile.open(out) as t:
        bedpe.write_bytes(t.extractfile("sample.bedpe.gz").read())
    records = sv.parse_bedpe(bedpe)
    assert records[0]["svclass"] == source_class  # the artifact preserves the source label
    records[0]["sample"] = "S1"
    g = pg.Grove(order=100)
    _, created = sv.attach_tracked(g, records)
    anchor = next(iter(g.intersect(pg.GenomicCoordinate("*", 0, 999), "chr1")))
    (_, edge), = g.get_edge_list(anchor)
    assert edge["svclass"] == expected_class
    assert edge["source_svclass"] == source_class
    assert edge["junction_class"] == expected_junction
    sv.detach(g, created)
    assert g.size() == 0 and g.edge_count() == 0
