# SPDX-License-Identifier: GPL-3.0-or-later
"""Derive the per-project PCAWG expression and copy-number tables from the pinned Xena
originals → two ready-to-pin tarballs and their catalogs.

    python tools/build_pcawg_xena.py <out_dir>

Build-time step, not run on a user's machine (same shape as `tools/build_hichip_cohorts.py`):

* `pcawg.rnaseq.xena` — the PCAWG-3 gene × specimen matrix (log2(FPKM-UQ + 0.001), GENCODE v19
  gene ids, 1,521 columns: ICGC specimen ids plus a few aliquot uuids). Re-cut into one wide
  table per ICGC project code, `<code>.tsv`: `gene` (ENSG, unversioned), `chrom`, `tss`
  (GENCODE v50, from the packaged `gene_tss.tsv.gz`), then one column per tumour named by its
  **WGS aliquot id** — the id the SV layer uses — resolved through the pinned sample sheet
  (specimen → the WGS tumour aliquot that is in the packaged SV catalog). Specimens with no
  such aliquot, and genes absent from the v50 backbone, are dropped and counted.
* `pcawg.cn.xena.hg19` — the PCAWG-11 consensus segments (specimen × segment, hg19). Lifted
  to GRCh38 with the same chain file `tools/liftover_pcawg_sv.py` uses, both ends
  independently; a segment is dropped if either end fails to lift, they land on different
  chromosomes, or the lifted end is not after the start. Re-cut per project with the same
  aliquot resolution, written BED half-open in the column shape `layers/cn.py` attaches.

Outputs: `pcawg_rnaseq_fpkm_uq_log2.by_cohort.tgz`, `pcawg_consensus_cn.hg38.by_cohort.tgz`,
`expression_cohorts.tsv`, `cn_cohorts.tsv` (the two catalogs to copy to
`src/genogrove_canopy/data/`), and each tarball's sha256. Requires `pyliftover`.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from genogrove_canopy import resources  # noqa: E402
from genogrove_canopy.layers.cn import _FIELDS as CN_FIELDS  # noqa: E402
from genogrove_canopy.layers.expression import _FIELDS as EXPR_FIELDS  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from liftover_pcawg_sv import _chain_file  # noqa: E402  the same pinned chain file


def _aliquot_map() -> tuple[dict, dict]:
    """``(specimen -> WGS aliquot, aliquot -> WGS aliquot)`` restricted to the aliquots in the
    packaged SV catalog, plus ``aliquot -> project code``."""
    sv = {a: r["project_code"] for r in resources.pcawg_cohorts() for a in r["aliquot_ids"]}
    with open(resources.resolve("pcawg.sample_sheet"), newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_specimen = {r["icgc_specimen_id"]: r["aliquot_id"] for r in rows if r["aliquot_id"] in sv}
    by_aliquot = {r["aliquot_id"]: by_specimen.get(r["icgc_specimen_id"]) for r in rows}
    return {**by_aliquot, **by_specimen}, sv


def _tss() -> dict:
    """``ensembl_base -> (chrom, tss)`` from the packaged GENCODE v50 table."""
    out = {}
    with gzip.open(Path(resources.__file__).with_name("data") / "gene_tss.tsv.gz", "rt") as fh:
        next(fh)
        for ln in fh:
            gid, _name, chrom, tss, _strand = ln.rstrip("\n").split("\t")
            out[gid] = (chrom, tss)
    return out


def build_expression(out_dir: Path, resolve: dict, project: dict) -> dict:
    tables = out_dir / "expression"
    tables.mkdir(parents=True, exist_ok=True)
    tss = _tss()
    by_code: dict[str, list[tuple[int, str]]] = defaultdict(list)  # code -> [(column index, aliquot)]
    unmatched = 0
    with open(resources.resolve("pcawg.rnaseq.xena")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for i, col in enumerate(header[1:], start=1):
            aliquot = resolve.get(col)
            if aliquot is None:
                unmatched += 1
                continue
            by_code[project[aliquot]].append((i, aliquot))
        handles = {}
        for code, cols in sorted(by_code.items()):
            handles[code] = open(tables / f"{code}.tsv", "w")
            handles[code].write("\t".join(EXPR_FIELDS + tuple(a for _, a in cols)) + "\n")
        n_genes, missed = 0, 0
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            gene = parts[0].split(".")[0]
            loc = tss.get(gene)
            if loc is None:
                missed += 1
                continue
            n_genes += 1
            for code, cols in by_code.items():
                handles[code].write("\t".join((gene, loc[0], loc[1]) + tuple(parts[i] for i, _ in cols)) + "\n")
    for h in handles.values():
        h.close()
    print(f"expression: {n_genes} genes kept, {missed} not in GENCODE v50; "
          f"{sum(len(c) for c in by_code.values())} tumours in {len(by_code)} cohorts, {unmatched} columns unmatched")
    return {code: (len(cols), n_genes) for code, cols in by_code.items()}


def build_cn(out_dir: Path, resolve: dict, project: dict, lo) -> dict:
    tables = out_dir / "cn"
    tables.mkdir(parents=True, exist_ok=True)
    handles, n_seg, samples, dropped, unmatched = {}, Counter(), defaultdict(set), Counter(), set()

    def lift(chrom, pos):
        hits = lo.convert_coordinate(chrom, pos)
        return hits[0][1] if hits and hits[0][0] == chrom and hits[0][2] == "+" else None

    with open(resources.resolve("pcawg.cn.xena.hg19"), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            aliquot = resolve.get(r["sampleID"])
            if aliquot is None:
                unmatched.add(r["sampleID"])
                continue
            code = project[aliquot]
            chrom = r["chr"] if r["chr"].startswith("chr") else f"chr{r['chr']}"
            start, end = lift(chrom, int(r["start"]) - 1), lift(chrom, int(r["end"]) - 1)  # 1-based inclusive -> 0-based
            if start is None or end is None or end < start:
                dropped[code] += 1
                continue
            out = handles.get(code)
            if out is None:
                out = handles[code] = open(tables / f"{code}.tsv", "w")
                out.write("\t".join(CN_FIELDS) + "\n")
            out.write("\t".join((chrom, str(start), str(end + 1), r["total_cn"], r["major_cn"],
                                 r["minor_cn"], r["star"], aliquot, code)) + "\n")
            n_seg[code] += 1
            samples[code].add(aliquot)
    for h in handles.values():
        h.close()
    print(f"cn: {sum(n_seg.values())} segments kept, {sum(dropped.values())} dropped in liftover, "
          f"{sum(len(s) for s in samples.values())} tumours in {len(handles)} cohorts, {len(unmatched)} specimens unmatched")
    return {code: (len(samples[code]), n_seg[code], dropped[code]) for code in handles}


def _tar(out_dir: Path, name: str, tables: Path) -> None:
    tgz = out_dir / name
    with tarfile.open(tgz, "w:gz") as dst:
        for table in sorted(tables.glob("*.tsv")):
            dst.add(table, arcname=table.name)
    print(f"{tgz}\t{tgz.stat().st_size} bytes\tsha256 {hashlib.sha256(tgz.read_bytes()).hexdigest()}")


def main() -> None:
    from pyliftover import LiftOver

    out_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xena_out")
    out_dir.mkdir(parents=True, exist_ok=True)
    resolve, project = _aliquot_map()

    expr = build_expression(out_dir, resolve, project)
    with (out_dir / "expression_cohorts.tsv").open("w") as fh:
        fh.write("project_code\tn_samples\tn_genes\n")
        for code, (n, g) in sorted(expr.items()):
            fh.write(f"{code}\t{n}\t{g}\n")
    _tar(out_dir, "pcawg_rnaseq_fpkm_uq_log2.by_cohort.tgz", out_dir / "expression")

    cn = build_cn(out_dir, resolve, project, LiftOver(str(_chain_file(out_dir))))
    with (out_dir / "cn_cohorts.tsv").open("w") as fh:
        fh.write("project_code\tn_samples\tn_segments\tn_dropped_liftover\n")
        for code, (n, s, d) in sorted(cn.items()):
            fh.write(f"{code}\t{n}\t{s}\t{d}\n")
    _tar(out_dir, "pcawg_consensus_cn.hg38.by_cohort.tgz", out_dir / "cn")


if __name__ == "__main__":
    main()
