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
# Derived, plain-text, one file per cohort — what the sandbox is granted and reads (it has no
# `gzip`). Rebuilt from the pinned index in under a second, so it is disposable cache, not an
# artifact: deleting it costs a rebuild, never correctness.
LINKS_DIR = resources._CACHE / "re2g_links"
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


def attach_links(grove, path, cohort):
    """Attach one cohort's links to a **mutable** ``grove``, from the plain table at ``path``.

    One **node per element** — an enhancer is one piece of DNA however many genes it regulates —
    and one ``regulates`` edge per (element, gene), carrying the evidence::

        node  {"type": "enhancer", "source": "ENCODE-rE2G", "class": "promoter|genic|intergenic"}
        edge  {"rel": "regulates", "byCohort": {<cohort>: {"score_max":.., "score_mean":.., "n_rep":..}}}

    The score belongs to the *link*, not to the interval, which is why it sits on the edge. The
    ``byCohort`` map is what makes a second cohort's call **merge** onto the same node and edge
    pair instead of adding a parallel node/edge: pygenogrove happily stores two edges between one
    pair, so the merge is done here — reuse the element node already at that interval, and fold
    the new cohort into the existing edge's map. ``regulated_by`` is stored as well so
    "its enhancers" from a gene is one plain ``get_neighbors_if`` hop.

    Returns ``(elements, links, missed)``.

    **This function is also shipped into the sandbox as source text** (see ``preamble``), so it
    must stay self-contained: everything it needs is imported inside it or passed in — no module
    globals, no helpers from this module. Keeping it a real function rather than a string literal
    is what lets ``tests/test_enhancers.py`` exercise the exact code the sandbox runs.
    """
    import pygenogrove as pg

    genes, nodes, links, missed = {}, {}, 0, 0
    with open(path) as fh:
        for ln in fh:
            chrom, start, end, cls, ens, tchrom, tss, n_rep, s_mean, s_max = \
                ln.rstrip("\n").split("\t")
            if ens not in genes:  # resolve the GENCODE gene once per gene, not per link
                hit = None
                # The table's TSS is 1-based (GFF); the grove is 0-based closed. Without the -1 a
                # minus-strand gene's TSS is its *end* plus one — one base outside the gene — and
                # every minus-strand target silently fails to resolve (110,403 of 225,229 links).
                for k in grove.intersect(pg.GenomicCoordinate("*", int(tss) - 1, int(tss) - 1),
                                         tchrom):
                    d = k.data
                    if d.get("type") == "gene" and (d.get("id") or "").split(".")[0] == ens:
                        hit = k
                        break
                genes[ens] = hit
            gene = genes[ens]
            if gene is None:  # target absent from this GENCODE build — count, never invent
                missed += 1
                continue
            ek = (chrom, int(start), int(end) - 1)  # rE2G BED half-open -> grove 0-based closed
            if ek not in nodes:  # one node per element — also across calls: a previous cohort
                nodes[ek] = next(  # may already have inserted this exact interval
                    (k for k in grove.intersect(pg.GenomicCoordinate("*", ek[1], ek[2]), chrom)
                     if k.data.get("source") == "ENCODE-rE2G"
                     and (k.value.start, k.value.end) == ek[1:]),
                    None,
                ) or grove.insert(chrom, pg.GenomicCoordinate(".", ek[1], ek[2]),
                                  {"type": "enhancer", "source": "ENCODE-rE2G", "class": cls})
            node = nodes[ek]
            # Edge payloads come back as copies, so a link already present from another cohort
            # is merged by removing the pair and re-adding it with the union of `byCohort`.
            by = next((m["byCohort"] for t, m in grove.get_edge_list(node)
                       if m and m.get("rel") == "regulates"
                       and t.data.get("id") == gene.data.get("id")), None)
            if by is not None:
                grove.remove_edge(node, gene)
                grove.remove_edge(gene, node)
            by = {**(by or {}), cohort: {"score_max": float(s_max), "score_mean": float(s_mean),
                                         "n_rep": int(n_rep)}}
            grove.add_edge(node, gene, {"rel": "regulates", "byCohort": by})
            grove.add_edge(gene, node, {"rel": "regulated_by", "byCohort": by})
            links += 1
    return len(nodes), links, missed


def links_file(cohort: str) -> Path:
    """The cohort's links as a plain TSV the sandbox can read, building it once if needed.

    Plain text, not the bgzip index itself, because the sandbox has no ``gzip`` — its allowlist is
    compute-only and its ``open`` is read-only and restricted to the roots the host grants. So the
    host decompresses the pinned index once per cohort (~0.6 s, ~18 MB) and grants that file.

    Columns are exactly what ``attach_links`` reads, in order: chrom, start, end, class, ensembl
    (unversioned), target chrom, target TSS, n_rep, score_mean, score_max. Rows whose target is
    absent from the gene table are dropped here rather than in the sandbox — resolving them there
    would fail anyway, and the host is where a count can be reported.
    """
    dest = LINKS_DIR / f"{_slug(cohort)}.links.tsv"
    if dest.exists():
        return dest
    if not ensure_index(cohort):
        raise KeyError(f"{cohort!r} is not in the pinned rE2G index bundle")
    by_ens, _ = _gene_tss()
    src = INDEX_DIR / f"{_slug(cohort)}.byEnhancer.tsv.gz"
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".tmp")
    with gzip.open(src, "rt") as fh, tmp.open("w") as out:
        for ln in fh:
            r = dict(zip(_FIELDS, ln.rstrip("\n").split("\t")))
            loc = by_ens.get(r["ensembl_id"].split(".")[0])
            if loc is None:
                continue
            out.write("\t".join((r["chrom"], r["start"], r["end"], r["class"],
                                 r["ensembl_id"].split(".")[0], loc[0], str(loc[1]),
                                 r["n_rep"], r["score_mean"], r["score_max"])) + "\n")
    tmp.replace(dest)  # atomic: a half-written table must never look cached
    return dest


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


def preamble(gg: str, cohort_links: dict[str, str] | None = None) -> str:
    """The sandbox preamble that binds ``GROVE`` to an open, ready-to-query handle.

    Without cohorts that is ``pg.GroveView.open(gg)`` — lazy, ~200 ms, exactly what a
    non-regulatory question paid before. With ``cohort_links`` (cohort id -> its links file path,
    see ``links_file``), the grove is deserialized **mutable** and each cohort's links are attached
    in the sandbox in turn (~6 s, ~1.8 GB for the first; an element linking to the same gene in two
    cohorts merges its ``byCohort`` payload onto one edge rather than duplicating — see
    ``attach_links``), because that is the only place the generated code can reach them: an
    in-memory grove cannot cross the process boundary, and the records cannot be injected as a
    literal (~34 MB of program text for one cohort alone).

    The build is memoised in ``_CANOPY_STATE``, which survives between queries in a warm worker,
    so an interactive session pays it once. The key is the grove path *and* the full set of
    (cohort, links) pairs, since a stale grove here is a silently wrong answer. Only one entry is
    kept: a different cohort selection evicts and rebuilds.

    ``ENHANCERS`` is still defined, and always empty. Generated code from an older prompt that
    loops over it gets nothing rather than a ``NameError``; the enhancers are in the grove now.
    """
    import json

    if not cohort_links:
        return f"import pygenogrove as pg\nGROVE = pg.GroveView.open({json.dumps(gg)})\nENHANCERS = []\n"

    import inspect

    attach_calls = "".join(
        f"    attach_links(GROVE, {json.dumps(links)}, {json.dumps(cohort)})\n"
        for cohort, links in cohort_links.items()
    )
    return (
        "import pygenogrove as pg\n"
        f"{inspect.getsource(attach_links)}\n"
        f"_key = ({json.dumps(gg)}, tuple(sorted({json.dumps(cohort_links)}.items())))\n"
        "_state = globals().get('_CANOPY_STATE')\n"       # absent in one-shot `sandbox.run`
        "if _state is not None and _state.get('key') == _key:\n"
        "    GROVE = _state['grove']\n"
        "else:\n"
        f"    GROVE = pg.Grove.deserialize({json.dumps(gg)})\n"
        f"{attach_calls}"
        "    if _state is not None:\n"
        "        _state.clear()\n"                        # one grove at a time; see the docstring
        "        _state.update(key=_key, grove=GROVE)\n"
        "ENHANCERS = []\n"
    )


def _index_files(cohort: str) -> list[str]:
    """The four filenames a cohort needs: two bgzip tables and their tabix indexes.

    Single source of truth for the set. When the presence check covered only the two `.tsv.gz`
    files, a cohort whose `.tbi` files had not arrived reported itself ready and then failed
    inside tabix.
    """
    slug = _slug(cohort)
    return [f"{slug}.{kind}{suffix}"
            for kind in ("byEnhancer.tsv.gz", "byTargetGene.tsv.gz")
            for suffix in ("", ".tbi")]


def index_present(cohort: str) -> bool:
    """True if **every** file the cohort needs is already cached — no fetching."""
    return all((INDEX_DIR / name).exists() for name in _index_files(cohort))


def ensure_index(cohort: str) -> bool:
    """Make ``cohort``'s index available, downloading the four pinned files if needed.

    Returns False when the cohort is not in the pinned bundle — a real answer ("we have no rE2G
    data for that biosample"), distinct from a failure. Anything else propagates: a checksum
    mismatch or a dead network should not quietly become "no enhancers found".
    """
    if index_present(cohort):
        return True
    names = _index_files(cohort)
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
         "whether a variant falls in an enhancer — scoped to a biosample/cohort",
    schema='enhancer node `{"type":"enhancer", "source":"ENCODE-rE2G", '
           '"class":<promoter|genic|intergenic>}` — filter on `source`; the target gene is one '
           '`get_neighbors` hop over the `{"rel":"regulates"}` edge (reverse: `regulated_by`), '
           'whose payload is `{"byCohort": {<cohort>: {"score_max":.., "score_mean":.., '
           '"n_rep":..}}}`. The score is on the **edge**, not the node — it is a property of the '
           'link, and one element may regulate several genes with different scores',
    attach=attach_links,
)
