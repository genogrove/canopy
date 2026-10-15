# SPDX-License-Identifier: GPL-3.0-or-later
"""Expression layer — one value per gene per tumour, hung off the gene node.

Source-agnostic: needs one wide table per cohort, a gene per row and a sample per column, the
cell being that sample's expression of that gene. The PCAWG-3 RNA-seq quantification (TopHat2
+ STAR averaged, FPKM-UQ, ``log2(x + 0.001)``, Xena ``tophat_star_fpkm_uq.v2_aliquot_gl.sp.log``)
is the validated example source, re-cut per PCAWG project by ``tools/build_pcawg_xena.py``.

The layer is a **payload**, not new intervals: a gene's expression is a property of the gene in
each tumour, so it is attached as one graph-only node per gene per cohort —

    node  {"type": "expression", "source": "PCAWG-RNAseq", "cohort": <project code>,
           "gene": <ENSG, unversioned>, "unit": "log2(fpkm-uq+0.001)", "values": {<sample>: v}}
    edge  gene -> node {"rel": "expressed",     "cohort": ..}
          node -> gene {"rel": "expression_of", "cohort": ..}

— ``add_external_key`` with the gene's own coordinate, so it is reachable only from the gene
and never comes back from ``intersect``. ``values`` is keyed by the tumour's **WGS aliquot id**,
the same id the SV layer puts on ``breakpoint_edge.sample``, so "carriers vs non-carriers" is a
plain dict lookup from the SV edges into this node; a specimen without a WGS aliquot in the SV
set is not in the table. One node per gene per cohort (not merged across cohorts) keeps the
attach to one ``insert`` per row and the filter on ``EXPRESSION_COHORTS`` a field test.

The table carries the gene's GENCODE v50 TSS (1-based, as ``gene_tss.tsv.gz`` has it) so the
sandbox can find the gene node by ``intersect`` at that base and match the unversioned id —
the same resolution the enhancer layer does. Genes absent from the backbone were dropped at
build time, never here.
"""

from __future__ import annotations

from pathlib import Path

from genogrove_canopy import resources
from genogrove_canopy.layers._base import Layer

#: The leading columns of a cohort table; every later column is one sample (WGS aliquot id).
_FIELDS = ("gene", "chrom", "tss")

EXPRESSION_DIR = resources._CACHE / "expression_cohorts"


def cohort_file(code: str) -> Path:
    """One PCAWG project's expression table as a plain TSV the sandbox can read, extracted once
    from the pinned derived tarball (``pcawg.expression``, one ``<code>.tsv`` per project, built
    by ``tools/build_pcawg_xena.py``) into ``<cache>/expression_cohorts/<pin>/<code>.tsv``."""
    import os
    import shutil
    import tarfile

    dest = (resources.pinned_dir(EXPRESSION_DIR, "pcawg.expression") / f"{code}.tsv").resolve()  # resolved: see cli._grove_context
    if dest.exists():
        return dest
    if code not in {r["project_code"] for r in resources.expression_cohorts()}:
        raise KeyError(f"{code!r} is not a PCAWG expression cohort (see --list-cohorts)")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.tmp")
    with tarfile.open(resources.resolve("pcawg.expression")) as t, tmp.open("wb") as out:
        shutil.copyfileobj(t.extractfile(f"{code}.tsv"), out)
    tmp.replace(dest)  # atomic: a half-written table must never look cached
    return dest


def attach(grove, path, cohort) -> int:
    """``Layer``-contract entry point: attach one cohort's table to a **mutable** grove and
    return the number of genes that got a node."""
    return attach_tracked(grove, path, cohort)[0]


def attach_tracked(grove, path, cohort):
    """Same as ``attach``, but also returns what it created, ``("node", gene_key, node_key)`` per
    gene, for ``detach``, plus the count of table genes not found in the backbone.

    **Shipped into the sandbox as source text** (see ``preamble.build``): self-contained —
    ``pygenogrove`` imported inside, no module globals, no helpers from this module.
    """
    import pygenogrove as pg

    created, n, missed = [], 0, 0
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[3:]
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            gene, chrom, tss = parts[0], parts[1], int(parts[2])
            # The table's TSS is 1-based (GFF); the grove is 0-based closed. Without the -1 a
            # minus-strand gene's TSS is its end plus one, one base outside the gene, and every
            # minus-strand gene silently fails to resolve (the enhancer layer's lesson).
            at = pg.GenomicCoordinate("*", tss - 1, tss - 1)
            hit = next((k for k in grove.intersect(at, chrom)
                        if k.data.get("type") == "gene" and (k.data.get("id") or "").split(".")[0] == gene), None)
            if hit is None:  # absent from this backbone — count, never invent
                missed += 1
                continue
            values = {s: float(v) for s, v in zip(samples, parts[3:]) if v not in ("", "NA")}
            node = grove.add_external_key(hit.value, {"type": "expression", "source": "PCAWG-RNAseq",
                                                      "cohort": cohort, "gene": gene,
                                                      "unit": "log2(fpkm-uq+0.001)", "values": values})
            grove.add_edge(hit, node, {"rel": "expressed", "cohort": cohort})
            grove.add_edge(node, hit, {"rel": "expression_of", "cohort": cohort})
            created.append(("node", hit, node))
            n += 1
    return n, created, missed


def detach(grove, created) -> None:
    """Remove this attachment's expression nodes and the two edges each carries. An external
    key cannot be removed from the pinned API, so the node is left edgeless and unreachable."""
    for _, gene, node in created:
        grove.remove_edge(gene, node)
        grove.remove_edge(node, gene)


LAYER = Layer(
    name="expression",
    axis="genomic",
    kind="payload",
    title="Gene expression — PCAWG RNA-seq per gene per tumour, per tumour cohort",
    when="a question asks whether a gene is expressed, up or down, or differently between "
         "tumours — carriers of a rearrangement vs the rest, one tumour vs its cohort — in one "
         "or more PCAWG tumour cohorts (project codes such as LIRI-JP); the sandbox filters on "
         "EXPRESSION_COHORTS",
    schema='one graph-only `{"type":"expression", "source":"PCAWG-RNAseq", "cohort":<PCAWG '
           'project code>, "gene":<ENSG>, "unit":"log2(fpkm-uq+0.001)", "values":{<WGS aliquot '
           'id>: <float>}}` node per gene per cohort, reached from the gene by `get_edge_list` '
           'filtered on `rel == "expressed"` (never returned by `intersect`). `values` is keyed '
           'by the same aliquot id the SV layer puts on `breakpoint_edge.sample`, so carriers '
           'vs non-carriers is a dict lookup: `{s: v for s, v in values.items() if s in '
           'carriers}` against the rest. Keep `cohort in EXPRESSION_COHORTS`. Report the unit, '
           'the two group sizes and a difference of medians — never a p-value you did not '
           'compute',
    attach=attach,
)
