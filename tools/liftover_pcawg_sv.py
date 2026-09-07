# SPDX-License-Identifier: GPL-3.0-or-later
"""Lift PCAWG's consensus SV BEDPE (hg19) to GRCh38 → two ready-to-pin tarballs.

Build-time input, not resolved on a user's machine — same shape as `encode.ccre.v4`:
the raw PCAWG tarballs (pinned as `pcawg.sv.icgc`/`pcawg.sv.tcga`, hg19) are the input;
this script's *output* (lifted to GRCh38, matching the pinned GENCODE v50 backbone) is
what actually gets re-hosted and pinned. No official hg38 PCAWG SV release exists, so
re-hosting a derived artifact is the only way to get one — unlike per-sample hosting
(tried and reverted elsewhere in this project's history), this is a real reason.

    python tools/liftover_pcawg_sv.py <out_dir>

Requires `pyliftover` (one-off tool dependency, not a package runtime dependency —
run `pip install pyliftover` before invoking this script) and the UCSC hg19->hg38
chain file (downloaded here, sha256-verified against a fixed pin).
"""

from __future__ import annotations

import gzip
import hashlib
import io
import sys
import tarfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from genogrove_canopy import resources  # noqa: E402

CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
CHAIN_SHA256 = "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"

_FIELDS = ("chrom1", "start1", "end1", "chrom2", "start2", "end2",
           "sv_id", "pe_support", "strand1", "strand2", "svclass", "svmethod")


def _chain_file(out_dir: Path) -> Path:
    dest = out_dir / "hg19ToHg38.over.chain.gz"
    if dest.exists():
        return dest
    urllib.request.urlretrieve(CHAIN_URL, dest)  # noqa: S310 — fixed UCSC host
    digest = hashlib.sha256(dest.read_bytes()).hexdigest()
    if digest != CHAIN_SHA256:
        dest.unlink()
        raise RuntimeError(f"chain file sha256 mismatch: got {digest}, expected {CHAIN_SHA256}")
    return dest


def _lift_one(lo, chrom: str, pos: int) -> int | None:
    """hg19 0-based ``pos`` -> hg38 0-based, or ``None`` if it doesn't lift cleanly
    (no result, or lands on a different/unplaced contig).

    PCAWG's own chrom column has no ``chr`` prefix (``"1"``, not ``"chr1"``); the
    chain file — and the GENCODE backbone this feeds into — both use ``"chr1"``.
    """
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    hits = lo.convert_coordinate(chrom, pos)
    if not hits or hits[0][0] != chrom:
        return None
    return hits[0][1]


def lift_tarball(src_tgz: Path, out_tgz: Path, lo) -> tuple[int, int]:
    """Lift every record in every sample file inside ``src_tgz``, writing a new
    tarball with the same structure at ``out_tgz``. Returns (kept, dropped)."""
    kept = dropped = 0
    with tarfile.open(src_tgz, "r:gz") as src, tarfile.open(out_tgz, "w:gz") as dst:
        for member in src.getmembers():
            if not member.name.endswith(".bedpe.gz"):
                continue
            fh = src.extractfile(member)
            with gzip.open(fh, "rt") as gz:
                header = gz.readline()
                lines = []
                for ln in gz:
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) != len(_FIELDS):
                        continue
                    row = dict(zip(_FIELDS, parts))
                    p1 = _lift_one(lo, row["chrom1"], int(row["start1"]))
                    p2 = _lift_one(lo, row["chrom2"], int(row["start2"]))
                    if p1 is None or p2 is None:
                        dropped += 1
                        continue
                    # Normalize to "chr1" form (PCAWG's own columns lack the prefix) so the
                    # output matches the GENCODE backbone's own chromosome naming.
                    row["chrom1"] = row["chrom1"] if row["chrom1"].startswith("chr") else f"chr{row['chrom1']}"
                    row["chrom2"] = row["chrom2"] if row["chrom2"].startswith("chr") else f"chr{row['chrom2']}"
                    row["start1"], row["end1"] = str(p1), str(p1 + 1)
                    row["start2"], row["end2"] = str(p2), str(p2 + 1)
                    lines.append("\t".join(row[f] for f in _FIELDS))
                    kept += 1
            out_bytes = gzip.compress((header + "\n".join(lines) + "\n").encode())
            info = tarfile.TarInfo(member.name)
            info.size = len(out_bytes)
            info.mtime = member.mtime
            dst.addfile(info, io.BytesIO(out_bytes))
    return kept, dropped


def main() -> None:
    from pyliftover import LiftOver

    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("liftover_out")
    out_dir.mkdir(parents=True, exist_ok=True)

    chain = _chain_file(out_dir)
    lo = LiftOver(str(chain))

    for name in ("pcawg.sv.icgc", "pcawg.sv.tcga"):
        src = resources.resolve(name)
        stem = Path(resources.RESOURCES[name].filename).stem  # strips the trailing .tgz
        out = out_dir / f"{stem}.hg38.tgz"
        kept, dropped = lift_tarball(src, out, lo)
        digest = hashlib.sha256(out.read_bytes()).hexdigest()
        print(f"{name}: {kept} SVs lifted, {dropped} dropped (failed to lift) "
              f"({100 * dropped / (kept + dropped):.2f}%)")
        print(f"  -> {out} ({out.stat().st_size / 1e6:.2f} MB, sha256={digest})")


if __name__ == "__main__":
    main()
