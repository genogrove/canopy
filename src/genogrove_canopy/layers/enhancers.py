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
from functools import lru_cache
from pathlib import Path

from genogrove_canopy import resources
from genogrove_canopy.layers._base import Layer

# The index bundle is pinned per file on Hugging Face and fetched a cohort at a time — four files,
# ~8 MB, not the whole 3 GB. Until it was pinned this directory had to be built locally, so a
# cohort absent from it made a lookup silently return nothing on every machine but the one that
# built it. gene_tss ships in the package (small).
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


def _slug(cohort: str) -> str:
    return cohort.replace(":", "_")


def attach_links(grove, path, cohort, nodes=None):
    """Attach one cohort's links to a **mutable** ``grove``, from the plain table at ``path``.

    One **node per element** — an enhancer is one piece of DNA however many genes it regulates —
    and one ``regulates`` edge per (element, gene), carrying the evidence::

        node  {"type": "enhancer", "source": "ENCODE-rE2G", "class": "promoter|genic|intergenic"}
        edge  {"rel": "regulates", "byCohort": {<cohort>: {"score_max":.., "score_mean":.., "n_rep":..}}}

    The score belongs to the *link*, not to the interval, which is why it sits on the edge. The
    ``byCohort`` map is what makes a second cohort's call **merge** onto the same node and edge
    pair instead of adding a parallel node/edge: pygenogrove happily stores two edges between one
    pair, so the merge is done here. ``nodes`` is the caller-owned element cache (interval ->
    Key) that makes it possible — pass the **same dict** for every cohort attached to one grove,
    as ``preamble`` does. A dict lookup is what keeps a cohort at ~8 s; finding the node by
    ``intersect`` instead cost 10 s more per cohort (408k lookups, measured). An element first
    seen in this call cannot yet carry another cohort's edge, so only shared elements pay the
    edge scan. ``regulated_by`` is stored as well so "its enhancers" from a gene is one plain
    ``get_neighbors_if`` hop.

    Returns ``(elements, links, missed)`` — ``elements`` is the size of ``nodes`` after the call.

    **This function is also shipped into the sandbox as source text** (see ``preamble``), so it
    must stay self-contained: everything it needs is imported inside it or passed in — no module
    globals, no helpers from this module. Keeping it a real function rather than a string literal
    is what lets ``tests/test_enhancers.py`` exercise the exact code the sandbox runs.
    """
    import pygenogrove as pg

    genes, links, missed = {}, 0, 0
    nodes = nodes if nodes is not None else {}
    new = set()  # elements this call inserted: no other cohort's edge can hang off them yet
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
            node, by = nodes.get(ek), None
            if node is None:  # one node per element, across every cohort attached to this grove
                node = nodes[ek] = grove.insert(chrom, pg.GenomicCoordinate(".", ek[1], ek[2]),
                                                {"type": "enhancer", "source": "ENCODE-rE2G",
                                                 "class": cls})
                new.add(ek)
            elif ek not in new:  # shared element: another cohort may already link it to this
                # gene (within one cohort a link is unique, so a node made here needs no scan).
                # payloads come back as copies, so merge = remove the pair, re-add the union.
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
    dest = (LINKS_DIR / f"{_slug(cohort)}.links.tsv").resolve()  # resolved: see cli._grove_context
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


def preamble(gg: str, cohort_links: dict[str, str] | None = None) -> str:
    """The sandbox preamble that binds ``GROVE`` to an open, ready-to-query handle.

    Without cohorts that is ``pg.GroveView.open(gg)`` — lazy, ~200 ms, exactly what a
    non-regulatory question paid before. With ``cohort_links`` (cohort id -> its links file path,
    see ``links_file``), the grove is deserialized **mutable** and each cohort's links are attached
    in the sandbox (~6 s to deserialize, ~8 s for a large cohort), because that is the only place
    the generated code can reach them: an in-memory grove cannot cross the process boundary, and
    the records cannot be injected as a literal (~34 MB of program text for one cohort alone).

    The grove is memoised in ``_CANOPY_STATE``, which survives between queries in a warm worker,
    **and grows**: a cohort asked about later is attached onto the same grove (the shared
    ``nodes`` cache is kept in the state so ``attach_links`` merges it onto existing elements and
    edges), never rebuilt. So after a breast-cancer question and then a prostate one, the grove
    holds both, one edge per shared link with both cohorts in its ``byCohort`` — which is what
    makes "in prostate but not breast" or "shared by both" a set operation on that map. The
    generated code is told which cohorts *this* question is about through ``COHORTS`` and must
    filter on it; the prompt says so.
    ponytail: the grove only grows within a session — a few hundred MB per cohort, no eviction.
    Add an LRU over ``cohorts`` if sessions that wander across many tissues turn up.
    """
    import json

    if not cohort_links:  # the worker hands _CANOPY_STATE to every query: hide it here too
        return (f"import pygenogrove as pg\nGROVE = pg.GroveView.open({json.dumps(gg)})\n"
                "COHORTS = []\nglobals().pop('_CANOPY_STATE', None)\n")

    import inspect

    return (
        "import pygenogrove as pg\n"
        f"{inspect.getsource(attach_links)}\n"
        "_state = globals().get('_CANOPY_STATE')\n"       # absent in one-shot `sandbox.run`
        "if _state is None:\n"
        "    _state = {}\n"
        f"if _state.get('gg') != {json.dumps(gg)}:\n"
        "    _state.clear()\n"   # drop a stale grove BEFORE deserializing the next one: both
        # live at once is ~1.8 GB + a ~2 GB deserialize peak against the 4 GiB sandbox cap
        f"    _state.update(gg={json.dumps(gg)}, grove=pg.Grove.deserialize({json.dumps(gg)}),\n"
        "                  nodes={}, cohorts=set())\n"
        "_grove = _state['grove']\n"
        f"for _c, _links in {json.dumps(sorted(cohort_links.items()))}:\n"
        "    if _c not in _state['cohorts']:\n"            # attach onto the warm grove, once
        "        attach_links(_grove, _links, _c, _state['nodes'])\n"
        "        _state['cohorts'].add(_c)\n"
        f"COHORTS = {json.dumps(sorted(cohort_links))}\n"  # this question's cohorts, by id
        # The memoised grove is a mutable `Grove` shared by every query of the session. Hand the
        # generated code a view that forwards only what a read-only `GroveView` has (plus the
        # count/lookup helpers), so a query cannot insert into it and quietly change the next
        # answer. A closure, not an attribute, so the raw grove isn't one `._g` away.
        "def _readonly(g):\n"
        "    _ok = {n for n in dir(pg.GroveView) if not n.startswith('_')} | {\n"
        "        'size', 'edge_count', 'vertex_count', 'external_vertex_count',\n"
        "        'indexed_vertex_count', 'vertex_count_with_edges', 'key_storage_size',\n"
        "        'has_edge', 'graph_empty'}\n"
        "    class _View:\n"
        "        def __getattr__(self, n):\n"
        "            if n in _ok:\n"
        "                return getattr(g, n)\n"
        "            raise AttributeError(f'GROVE is read-only: {n!r} is not a query method')\n"
        "        def __len__(self):\n"
        "            return len(g)\n"
        "    return _View()\n"
        "GROVE = _readonly(_grove)\n"
        # Drop the host-only names from the namespace the generated code runs in, so a query
        # won't by accident evict or swap the grove the next query in a warm session is given.
        # (`sys.modules['__main__']` can still reach it — the sandbox's documented residual risk.)
        "for _n in ('_CANOPY_STATE', '_state', '_grove', '_c', '_links', '_readonly', '_n'):\n"
        "    globals().pop(_n, None)\n"
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
