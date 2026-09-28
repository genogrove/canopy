# SPDX-License-Identifier: GPL-3.0-or-later
"""Derive the per-project TCGA HiChIP loop tables from the pinned GDC originals → one
ready-to-pin tarball.

    python tools/build_hichip_cohorts.py <out_dir>

Build-time step, not run on a user's machine (same shape as `tools/liftover_pcawg_sv.py`):
the pinned `tcga.hichip.metadata` manifest maps each FitHiChIP file to its TCGA project and
aliquot; the pinned `tcga.hichip.loops.raw` tarball holds one loop file per sample. This
script writes `<out_dir>/tcga_hichip_fithichip_loops_Q0.1.hg38.tgz` containing one
`<project_id>.tsv` per TCGA project (`TCGA-BRCA.tsv`, …) in the column shape
`layers/hichip.py` attaches — the layer's `cohort_file` extracts one project's table from it
at query time, as `layers/sv.py` does from the PCAWG tarballs — plus `hichip_cohorts.tsv`,
the per-project catalog to copy to `src/genogrove_canopy/data/`.

No re-thresholding and no lifting: the calls stay exactly as the source made them (q <= 0.1,
already GRCh38); only the 19 per-bin QC columns FitHiChIP carries are dropped. Bins stay
half-open BED as in the source — the layer converts to the grove's 0-based closed form.
"""

from __future__ import annotations

import csv
import hashlib
import sys
import tarfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from genogrove_canopy import resources  # noqa: E402
from genogrove_canopy.layers.hichip import _FIELDS  # noqa: E402

# FitHiChIP column -> derived column. Everything else in the source row is per-bin QC.
_SOURCE = {"chr1": "chrom1", "s1": "start1", "e1": "end1", "chr2": "chrom2", "s2": "start2",
           "e2": "end2", "cc": "cc", "Q-Value_Bias": "q", "isPeak1": "peak1", "isPeak2": "peak2"}


def main() -> None:
    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("hichip_out")
    tables = out_dir / "tcga-hichip"
    tables.mkdir(parents=True, exist_ok=True)

    with open(resources.resolve("tcga.hichip.metadata"), newline="") as fh:
        sample = {r["file_name"]: (r["project_id"], r["aliquot_submitter_id"])
                  for r in csv.DictReader(fh, delimiter="\t")
                  if r["tar_archive"] == "TCGA_HiChIP_FitHiChIP_loop_calls.tar.gz"}

    n_loops, n_samples, handles = Counter(), Counter(), {}
    with tarfile.open(resources.resolve("tcga.hichip.loops.raw"), "r:gz") as t:
        for m in t.getmembers():
            if not m.isfile():
                continue
            project, aliquot = sample[m.name.split("/")[-1]]  # KeyError = a file the manifest lacks: stop
            out = handles.get(project)
            if out is None:
                out = handles[project] = open(tables / f"{project}.tsv", "w")
                out.write("\t".join(_FIELDS) + "\n")
            n_samples[project] += 1
            fh = t.extractfile(m)
            header = fh.readline().decode().rstrip("\n").split("\t")
            idx = [header.index(src) for src in _SOURCE]  # the manifest guarantees one layout
            for ln in fh:
                parts = ln.decode().rstrip("\n").split("\t")
                row = dict(zip(_SOURCE.values(), (parts[i] for i in idx)))
                out.write("\t".join([row[f] for f in _FIELDS[:-2]] + [aliquot, project]) + "\n")
                n_loops[project] += 1
    for out in handles.values():
        out.close()

    tgz = out_dir / "tcga_hichip_fithichip_loops_Q0.1.hg38.tgz"
    with tarfile.open(tgz, "w:gz") as dst:
        for project in sorted(handles):
            dst.add(tables / f"{project}.tsv", arcname=f"{project}.tsv")
    # The packaged catalog (`src/genogrove_canopy/data/hichip_cohorts.tsv`): what
    # `--list-cohorts` shows and `cohort_file` validates against, without the tarball.
    with (out_dir / "hichip_cohorts.tsv").open("w") as fh:
        fh.write("project_code\tn_samples\tn_loops\n")
        for project in sorted(handles):
            fh.write(f"{project}\t{n_samples[project]}\t{n_loops[project]}\n")
    for project in sorted(handles):
        print(f"{project}\t{n_samples[project]} samples\t{n_loops[project]} loops")
    print(f"{tgz}\t{tgz.stat().st_size} bytes\tsha256 {hashlib.sha256(tgz.read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
