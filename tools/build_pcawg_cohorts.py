# SPDX-License-Identifier: GPL-3.0-or-later
"""Derive the packaged PCAWG cohort catalog from the pinned sample sheet + SV tarballs.

    python tools/build_pcawg_cohorts.py > src/genogrove_canopy/data/pcawg_cohorts.tsv

One row per ICGC project code that has at least one SV sample: code, number of SV samples,
number of SV records, and the aliquot ids (comma-separated) so the host can extract a cohort
from the tarballs without reading the sample sheet at query time. Build-time only: every
input is a pinned Resource, so the output is reproducible.
"""

from __future__ import annotations

import csv
import gzip
import io
import sys
import tarfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from genogrove_canopy import resources  # noqa: E402


def main() -> None:
    with open(resources.resolve("pcawg.sample_sheet"), newline="") as fh:
        project = {r["aliquot_id"]: r["dcc_project_code"] for r in csv.DictReader(fh, delimiter="\t")}
    members, n_sv = defaultdict(list), defaultdict(int)
    for name in ("pcawg.sv.icgc", "pcawg.sv.tcga"):
        with tarfile.open(resources.resolve(name)) as t:
            for m in t.getmembers():
                if not m.name.endswith(".bedpe.gz"):
                    continue
                aliquot = m.name.split("/")[-1].split(".")[0]
                code = project[aliquot]  # KeyError = a sample the sheet does not know: stop
                members[code].append(aliquot)
                n_sv[code] += sum(1 for _ in gzip.open(io.BytesIO(t.extractfile(m).read()), "rt")) - 1
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(["project_code", "n_samples", "n_svs", "aliquot_ids"])
    for code in sorted(members):
        out.writerow([code, len(members[code]), n_sv[code], ",".join(sorted(members[code]))])


if __name__ == "__main__":
    main()
