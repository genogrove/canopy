# SPDX-License-Identifier: GPL-3.0-or-later
"""Derive the per-project PCAWG expression and copy-number tables from the pinned Xena
originals → two ready-to-pin tarballs and their catalogs.

    python tools/build_pcawg_xena.py <out_dir>

Build-time step, not run on a user's machine (same shape as `tools/build_hichip_cohorts.py`):

* `pcawg.rnaseq.xena` — the PCAWG-3 gene × specimen matrix (log2(FPKM-UQ + 0.001), GENCODE v19
  gene ids, 1,521 columns: 1,466 ICGC specimen ids plus 55 aliquot uuids the sample sheet does
  not know, which stay unmatched). Re-cut into one wide table per ICGC project code,
  `<code>.tsv`: `gene` (ENSG, unversioned), `chrom`, `tss` (GENCODE v50, 1-based, from the
  packaged `gene_tss.tsv.gz`), then one column per tumour named by its **WGS aliquot id** —
  the id the SV layer uses — resolved through the pinned sample sheet: the specimen's own WGS
  tumour aliquot when it is in the packaged SV catalog, else — RNA-seq is often registered
  under a different specimen id than the WGS of the same tumour — the donor's SV aliquot when
  the donor has exactly one. Specimens that resolve neither way, and genes absent from the
  v50 backbone, are dropped and counted.
* `pcawg.cn.xena.hg19` — the PCAWG-11 consensus segments (specimen × segment, hg19). Lifted
  to GRCh38 with the same chain file `tools/liftover_pcawg_sv.py` uses, both ends
  independently; a segment is dropped if either end fails to lift, lands on another
  chromosome or in a reversed block, or the lifted end is not after the start
  (`n_dropped_liftover`), and also if the lifted length differs from the original by more
  than 10 % **and** more than 1 Mb (`n_dropped_length`) — two ends that each lift cleanly can
  still bracket inserted or removed sequence, and a stretched segment would assert a copy
  number over bases the call never covered. Re-cut per project with the same aliquot
  resolution, written BED half-open in the column shape `layers/cn.py` attaches.

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


def resolve_specimens(rows: list[dict], sv: dict) -> tuple[dict, set]:
    """``(specimen id -> WGS aliquot id, specimens that are the aliquot's own)`` for the aliquots
    in ``sv`` (aliquot -> project code), from sample-sheet ``rows``: a specimen's own SV aliquot
    first; else, for a tumour specimen of a donor with exactly one SV aliquot, that aliquot
    (RNA-seq and WGS of one tumour are often registered under different specimen ids). A donor
    with several SV tumours gets no guess — the RNA could belong to any of them. The caller
    uses the fallback only for an aliquot that has no own column (see ``build_expression``)."""
    own = {r["icgc_specimen_id"]: r["aliquot_id"] for r in rows if r["aliquot_id"] in sv}
    donor_aliquots = defaultdict(set)
    for r in rows:
        if r["aliquot_id"] in sv:
            donor_aliquots[r["icgc_donor_id"]].add(r["aliquot_id"])
    out = dict(own)
    for r in rows:
        spec = r["icgc_specimen_id"]
        if spec in out or r["dcc_specimen_type"].startswith("Normal"):
            continue
        aliquots = donor_aliquots.get(r["icgc_donor_id"], ())
        if len(aliquots) == 1:
            out[spec] = next(iter(aliquots))
    return out, set(own)


def _aliquot_map() -> tuple[dict, dict]:
    sv = {a: r["project_code"] for r in resources.pcawg_cohorts() for a in r["aliquot_ids"]}
    with open(resources.resolve("pcawg.sample_sheet"), newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    mapping, own = resolve_specimens(rows, sv)
    return mapping, own, sv


#: A lifted segment may differ from the original by this much before it is dropped: more
#: than both the fraction and the absolute length. Inserted hg38 centromere models stretch a
#: few long segments legitimately; a 10 Mb call becoming 33 Mb is not that.
LENGTH_TOLERANCE = (0.10, 1_000_000)


def lift_segment(lo, chrom: str, start: int, end: int):
    """hg19 1-based inclusive ``[start, end]`` -> hg38 0-based half-open ``(start, end)``, or
    ``None`` (``"liftover"``: an end fails, lands elsewhere, in a reversed block, or before the
    start) or ``"length"`` (both lift, but the segment stretched or shrank beyond
    :data:`LENGTH_TOLERANCE`). Chain positions are 0-based; a reversed block is refused because
    the segment's extent on the other strand is not what its two ends say."""
    def one(pos):
        hits = lo.convert_coordinate(chrom, pos)
        return hits[0][1] if hits and hits[0][0] == chrom and hits[0][2] == "+" else None

    s, e = one(start - 1), one(end - 1)
    if s is None or e is None or e < s:
        return None
    before, after = end - start + 1, e - s + 1
    frac, absolute = LENGTH_TOLERANCE
    if abs(after - before) > frac * before and abs(after - before) > absolute:
        return "length"
    return s, e + 1


def _tss() -> dict:
    """``ensembl_base -> (chrom, tss)`` from the packaged GENCODE v50 table."""
    out = {}
    with gzip.open(Path(resources.__file__).with_name("data") / "gene_tss.tsv.gz", "rt") as fh:
        next(fh)
        for ln in fh:
            gid, _name, chrom, tss, _strand = ln.rstrip("\n").split("\t")
            out[gid] = (chrom, tss)
    return out


def pick_columns(header: list[str], resolve: dict, own: set) -> tuple[dict, int, int]:
    """``aliquot -> column index`` for the matrix ``header`` (first cell is the feature name):
    the aliquot's own specimen column when there is one; else the donor-fallback column, and
    only if it is the single such column — a donor with a primary and two metastases has three
    RNA columns for one SV aliquot, and none of them is known to be the sequenced tumour.
    Returns the picks, the unmatched column count and the dropped-duplicate count."""
    chosen, fallback, unmatched = {}, defaultdict(list), 0
    for i, col in enumerate(header[1:], start=1):
        aliquot = resolve.get(col)
        if aliquot is None:
            unmatched += 1
        elif col in own:
            chosen[aliquot] = i
        else:
            fallback[aliquot].append(i)
    dropped = 0
    for aliquot, cols in fallback.items():
        if aliquot in chosen or len(cols) > 1:
            dropped += len(cols)
        else:
            chosen[aliquot] = cols[0]
    return chosen, unmatched, dropped


def build_expression(out_dir: Path, resolve: dict, own: set, project: dict) -> dict:
    tables = out_dir / "expression"
    tables.mkdir(parents=True, exist_ok=True)
    tss = _tss()
    by_code: dict[str, list[tuple[int, str]]] = defaultdict(list)  # code -> [(column index, aliquot)]
    with open(resources.resolve("pcawg.rnaseq.xena")) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        chosen, unmatched, duplicates = pick_columns(header, resolve, own)
        for aliquot, i in sorted(chosen.items(), key=lambda kv: kv[1]):
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
          f"{sum(len(c) for c in by_code.values())} tumours in {len(by_code)} cohorts, {unmatched} columns unmatched, "
          f"{duplicates} duplicate columns dropped")
    return {code: (len(cols), n_genes) for code, cols in by_code.items()}


def build_cn(out_dir: Path, resolve: dict, project: dict, lo) -> dict:
    tables = out_dir / "cn"
    tables.mkdir(parents=True, exist_ok=True)
    handles, n_seg, samples, unmatched = {}, Counter(), defaultdict(set), set()
    dropped = {"liftover": Counter(), "length": Counter()}
    with open(resources.resolve("pcawg.cn.xena.hg19"), newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            aliquot = resolve.get(r["sampleID"])
            if aliquot is None:
                unmatched.add(r["sampleID"])
                continue
            code = project[aliquot]
            chrom = r["chr"] if r["chr"].startswith("chr") else f"chr{r['chr']}"
            lifted = lift_segment(lo, chrom, int(r["start"]), int(r["end"]))
            if lifted is None or lifted == "length":
                dropped["liftover" if lifted is None else "length"][code] += 1
                continue
            start, end = lifted
            out = handles.get(code)
            if out is None:
                out = handles[code] = open(tables / f"{code}.tsv", "w")
                out.write("\t".join(CN_FIELDS) + "\n")
            out.write("\t".join((chrom, str(start), str(end), r["total_cn"], r["major_cn"],
                                 r["minor_cn"], r["star"], aliquot, code)) + "\n")
            n_seg[code] += 1
            samples[code].add(aliquot)
    for h in handles.values():
        h.close()
    print(f"cn: {sum(n_seg.values())} segments kept, {sum(dropped['liftover'].values())} dropped in liftover, "
          f"{sum(dropped['length'].values())} dropped for length, "
          f"{sum(len(s) for s in samples.values())} tumours in {len(handles)} cohorts, {len(unmatched)} specimens unmatched")
    return {code: (len(samples[code]), n_seg[code], dropped["liftover"][code], dropped["length"][code]) for code in handles}


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
    resolve, own, project = _aliquot_map()

    expr = build_expression(out_dir, resolve, own, project)
    with (out_dir / "expression_cohorts.tsv").open("w") as fh:
        fh.write("project_code\tn_samples\tn_genes\n")
        for code, (n, g) in sorted(expr.items()):
            fh.write(f"{code}\t{n}\t{g}\n")
    _tar(out_dir, "pcawg_rnaseq_fpkm_uq_log2.by_cohort.tgz", out_dir / "expression")

    cn = build_cn(out_dir, resolve, project, LiftOver(str(_chain_file(out_dir))))
    with (out_dir / "cn_cohorts.tsv").open("w") as fh:
        fh.write("project_code\tn_samples\tn_segments\tn_dropped_liftover\tn_dropped_length\n")
        for code, (n, s, d, dl) in sorted(cn.items()):
            fh.write(f"{code}\t{n}\t{s}\t{d}\t{dl}\n")
    _tar(out_dir, "pcawg_consensus_cn.hg38.by_cohort.tgz", out_dir / "cn")


if __name__ == "__main__":
    main()
