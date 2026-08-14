# SPDX-License-Identifier: GPL-3.0-or-later
"""Enhancer access layer — the fast, targeted replacement for a full per-cohort grove augment.

The sandbox that runs generated query code has a strict import allowlist (``pygenogrove`` +
compute-only modules — no ``os``/``subprocess``/``gzip``/file I/O), so the query code can only
see enhancers **through the grove**. This module is therefore **host-side**: it uses two tabix
indexes per cohort to fetch *only the enhancers a question needs*, which the host then attaches
to the GENCODE grove (``attach_to_grove``) before running the query — instead of materialising
a cohort's ~100k enhancers into a grove up front.

Indexes (built by tools, one pair per cohort, in ``INDEX_DIR``):
  <cohort>.byEnhancer.tsv.gz    — keyed by enhancer coordinate  -> region / variant queries
  <cohort>.byTargetGene.tsv.gz  — keyed by target-gene TSS       -> "enhancers of <gene>"
``gene_tss.tsv.gz`` (bundled in the package) maps gene symbol/ENSG -> TSS, parsed in ~65 ms.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

from genogrove_canopy import resources
from genogrove_canopy.layers._base import Layer

# The index bundle is pinned per file on Hugging Face and fetched a cohort at a time — four files,
# ~8 MB, not the whole 3 GB. Until it was pinned this directory had to be built locally, so a
# cohort absent from it made `fetch_for_targets` silently return nothing on every machine but the
# one that built it. gene_tss ships in the package (small).
INDEX_DIR = resources._CACHE / "re2g_index"
# Package data lives in genogrove_canopy/data/; this module is in genogrove_canopy/layers/,
# so go up one level.
_GENE_TSS = Path(__file__).parent.parent / "data" / "gene_tss.tsv.gz"

# Per-row columns of the consolidated edge (the enhancer→gene record the indexes carry).
_FIELDS = ("chrom", "start", "end", "target_gene", "ensembl_id", "class",
           "is_self_promoter", "cohort", "n_rep", "score_mean", "score_max")


@lru_cache(maxsize=1)
def _gene_tss():
    """``({ensembl_base: (chrom, tss)}, {symbol: ensembl_base})`` from the bundled TSV (~65 ms)."""
    by_ens, by_sym = {}, {}
    with gzip.open(_GENE_TSS, "rt") as fh:
        next(fh)  # header
        for ln in fh:
            gid, name, chrom, tss, _strand = ln.rstrip("\n").split("\t")
            by_ens[gid] = (chrom, int(tss))
            if name:
                by_sym[name] = gid
    return by_ens, by_sym


def resolve_gene(gene: str):
    """A gene symbol or ENSG id -> ``(ensembl_base, chrom, tss)``, or ``None`` if unknown."""
    by_ens, by_sym = _gene_tss()
    ens = (gene if gene.upper().startswith("ENSG") else by_sym.get(gene) or "").split(".")[0]
    loc = by_ens.get(ens)
    return (ens, *loc) if loc else None


def _slug(cohort: str) -> str:
    return cohort.replace(":", "_")


def _require_tabix() -> None:
    if shutil.which("tabix") is None:
        raise RuntimeError("`tabix` not found — install htslib (e.g. `brew install htslib`).")


def _query(index: Path, region: str, key_cols: int) -> list[dict]:
    """tabix ``region`` on ``index``; drop ``key_cols`` leading key columns; rows -> dicts."""
    if not index.exists():
        raise FileNotFoundError(f"{index} missing — enhancer index not present for this cohort "
                                f"(expected under {INDEX_DIR}).")
    _require_tabix()
    out = subprocess.run(["tabix", str(index), region],
                         capture_output=True, text=True, check=True).stdout
    recs = []
    for ln in out.splitlines():
        if not ln:
            continue
        f = ln.split("\t")[key_cols:]
        recs.append(dict(zip(_FIELDS, f)))
    return recs


def enhancers_of_gene(gene: str, cohort: str) -> list[dict]:
    """Enhancers predicted to regulate ``gene`` in ``cohort`` (Index B, by target-gene TSS).

    ``gene`` is a symbol (``"MYC"``) or ENSG id; ``cohort`` is a biosample ontology id
    (``"EFO:0005726"``). Empty list if the gene is unknown or has no linked enhancers.
    """
    hit = resolve_gene(gene)
    if hit is None:
        return []
    _ens, chrom, tss = hit
    return _query(INDEX_DIR / f"{_slug(cohort)}.byTargetGene.tsv.gz", f"{chrom}:{tss}-{tss}", 3)


def enhancers_in_region(chrom: str, start: int, end: int, cohort: str) -> list[dict]:
    """Enhancers overlapping ``chrom:start-end`` in ``cohort`` (Index A, by enhancer coordinate),
    each carrying its ``target_gene`` — for ``variant ∩ enhancer -> gene`` questions."""
    return _query(INDEX_DIR / f"{_slug(cohort)}.byEnhancer.tsv.gz", f"{chrom}:{start}-{end}", 0)


def attach_to_grove(grove, records, attached=None) -> int:
    """Augment a **mutable** GENCODE ``grove`` with ``records`` as first-class enhancer nodes +
    a ``regulates``/``regulated_by`` edge pair to the target gene node — so the query code
    traverses variant → enhancer → gene → (its exons / other enhancers) natively.

    Enhancers are **spatially inserted** (``g.insert``, not ``add_external_key``) so a variant
    ``intersect`` finds them alongside genes. The target gene node is located by ``intersect`` at
    the enriched ``target_chrom``/``target_tss`` and matched on Ensembl id. ``attached`` is a
    caller-owned set for **incremental** augmentation across a warm session — records already in
    it are skipped (no reset needed); returns the number newly attached.

    This is the targeted replacement for ``resources.augment_grove``: only the enhancers a query
    needs, onto an already-warm mutable grove — milliseconds, not a whole-cohort rebuild.
    """
    import pygenogrove as pg

    attached = attached if attached is not None else set()
    gene_cache: dict = {}
    n = 0
    for r in records:
        key = (r["cohort"], r["chrom"], r["start"], r["end"], r["target_gene"])
        if key in attached:
            continue
        ens = r["ensembl_id"].split(".")[0]
        gene = gene_cache.get(ens, False)
        if gene is False:  # resolve the GENCODE gene node once per gene
            gene = None
            tc, tt = r.get("target_chrom"), r.get("target_tss")
            if tc and tt:
                for k in grove.intersect(pg.GenomicCoordinate("*", int(tt), int(tt)), tc):
                    d = k.data
                    if d.get("type") == "gene" and (d.get("id") or "").split(".")[0] == ens:
                        gene = k
                        break
            gene_cache[ens] = gene
        # BED half-open -> grove 0-based closed (end - 1); unstranded enhancer.
        enh = grove.insert(r["chrom"], pg.GenomicCoordinate(".", int(r["start"]), int(r["end"]) - 1),
                           {"type": "enhancer", "class": r["class"], "target": r["target_gene"],
                            "score": float(r["score_max"]), "n_rep": int(r["n_rep"]),
                            "cohort": r["cohort"]})
        if gene is not None:
            meta = {"cohort": r["cohort"], "score": float(r["score_max"]), "n_rep": int(r["n_rep"])}
            grove.add_edge(enh, gene, {"rel": "regulates", **meta})
            grove.add_edge(gene, enh, {"rel": "regulated_by", **meta})
        attached.add(key)
        n += 1
    return n


def fetch_for_targets(targets, cohorts) -> list[dict]:
    """Enhancers for the LLM-declared ``targets`` across the selected ``cohorts`` (ontology ids).

    ``targets`` is a list of ``{"gene": "MYC"}`` and/or ``{"region": "chr8:127000000-128000000"}``.
    Deduped by (coord, target gene, cohort). This is what the host injects as ``ENHANCERS`` before
    running the generated query — only the enhancers the question actually asked for.
    """
    by_ens, _ = _gene_tss()
    out, seen = [], set()
    for cohort in cohorts:
        if not ensure_index(cohort):
            continue
        for t in targets or []:
            if t.get("gene"):
                recs = enhancers_of_gene(t["gene"], cohort)
            elif t.get("region"):
                chrom, span = t["region"].split(":")
                s, e = span.replace(",", "").split("-")
                recs = enhancers_in_region(chrom, int(s), int(e), cohort)
            else:
                continue
            for r in recs:
                key = (r["chrom"], r["start"], r["end"], r["target_gene"], r["cohort"])
                if key in seen:
                    continue
                seen.add(key)
                # Enrich with the target gene's locus so a grove augment can place the
                # regulates/regulated_by edge (and a query can hop into the gene's structure).
                loc = by_ens.get(r["ensembl_id"].split(".")[0])
                if loc:
                    r["target_chrom"], r["target_tss"] = loc[0], str(loc[1])
                out.append(r)
    return out


def preamble(records) -> str:
    """The sandbox preamble line that injects ``records`` as the ``ENHANCERS`` variable.

    Emitted as ``json.loads(<literal>)`` (``json`` is on the sandbox allowlist) so booleans/nulls
    round-trip safely rather than relying on JSON literals being valid Python."""
    import json

    return f"ENHANCERS = json.loads({json.dumps(json.dumps(records))})\n"


def index_present(cohort: str) -> bool:
    """True if both index files are already cached for ``cohort`` — no fetching."""
    s = _slug(cohort)
    return (INDEX_DIR / f"{s}.byTargetGene.tsv.gz").exists() and \
           (INDEX_DIR / f"{s}.byEnhancer.tsv.gz").exists()


def ensure_index(cohort: str) -> bool:
    """Make ``cohort``'s index available, downloading the four pinned files if needed.

    Returns False when the cohort is not in the pinned bundle — a real answer ("we have no rE2G
    data for that biosample"), distinct from a failure. Anything else propagates: a checksum
    mismatch or a dead network should not quietly become "no enhancers found".
    """
    if index_present(cohort):
        return True
    slug = _slug(cohort)
    names = [f"{slug}.{kind}{suffix}"
             for kind in ("byEnhancer.tsv.gz", "byTargetGene.tsv.gz")
             for suffix in ("", ".tbi")]
    manifest = resources.re2g_index_manifest()
    if any(n not in manifest for n in names):
        return False
    for name in names:
        resources.re2g_index_file(name)
    return True


LAYER = Layer(
    name="enhancers",
    axis="genomic",
    kind="edge",
    title="ENCODE-rE2G enhancer→gene links",
    when="a question is about enhancers, gene regulation, which enhancers regulate a gene, or "
         "whether a variant falls in an enhancer — optionally scoped to a biosample/cohort",
    schema='enhancer node `{"type":"enhancer", "class":.., "target":.., "score":.., "cohort":..}` '
           'plus a `{"rel":"regulates"}` / `{"rel":"regulated_by"}` edge pair to the GENCODE gene '
           '(each carrying `cohort`/`score`/`n_rep`)',
    attach=attach_to_grove,
)
