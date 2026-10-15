# SPDX-License-Identifier: GPL-3.0-or-later
"""Copy-number layer — consensus segments as indexed intervals, per tumour.

Source-agnostic: needs segments as half-open BED intervals with a total and, when the caller
gives them, major/minor allele copy numbers, per sample. The PCAWG-11 consensus copy number
(six callers reconciled on the consensus SV breakpoints; Xena ``20170119_final_consensus_
copynumber_sp``, hg19) is the validated example source, lifted to GRCh38 and re-cut per PCAWG
project by ``tools/build_pcawg_xena.py``.

A segment **is** an interval, so it is one indexed key per segment per sample —

    key  {"type": "cn_segment", "source": "PCAWG-consensus-CN", "sample": <WGS aliquot id>,
          "cohort": <project code>, "total_cn": int, "major_cn": int|null, "minor_cn": int|null,
          "star": int|null}

— returned by ``intersect`` alongside genes and cCREs, so "the copy number of gene G in sample
S" is ``intersect`` on the gene span filtered on type and sample. No edges: a segment relates
to a gene by overlap only. ``sample`` is the tumour's WGS aliquot id, the same id the SV layer
puts on ``breakpoint_edge.sample`` and the expression layer keys ``values`` by, so the three
join on it. Inserted additively per cohort onto the warm grove; ``detach`` removes exactly the
keys an attach made.

Segment ends were lifted independently; a segment whose ends land on different chromosomes,
reverse, or do not lift was dropped at build time (counted in the catalog), never guessed.
"""

from __future__ import annotations

from pathlib import Path

from genogrove_canopy import resources
from genogrove_canopy.layers._base import Layer

# One record per segment per sample, BED half-open, as the derived cohort table writes them.
_FIELDS = ("chrom", "start", "end", "total_cn", "major_cn", "minor_cn", "star", "sample", "cohort")

CN_DIR = resources._CACHE / "cn_cohorts"


def cohort_file(code: str) -> Path:
    """One PCAWG project's segments as a plain TSV the sandbox can read, extracted once from
    the pinned derived tarball (``pcawg.cn``, one ``<code>.tsv`` per project, built by
    ``tools/build_pcawg_xena.py``) into ``<cache>/cn_cohorts/<code>.tsv``."""
    import os
    import shutil
    import tarfile

    dest = (CN_DIR / f"{code}.tsv").resolve()  # resolved: see cli._grove_context
    if dest.exists():
        return dest
    if code not in {r["project_code"] for r in resources.cn_cohorts()}:
        raise KeyError(f"{code!r} is not a PCAWG copy-number cohort (see --list-cohorts)")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.tmp")
    with tarfile.open(resources.resolve("pcawg.cn")) as t, tmp.open("wb") as out:
        shutil.copyfileobj(t.extractfile(f"{code}.tsv"), out)
    tmp.replace(dest)  # atomic: a half-written table must never look cached
    return dest


def attach(grove, records) -> int:
    """``Layer``-contract entry point: attach segments to a **mutable** grove, return their count."""
    return attach_tracked(grove, records)[0]


def attach_tracked(grove, records):
    """Same as ``attach``, but also returns ``("key", chrom, Key)`` per segment inserted, for
    ``detach``.

    **Shipped into the sandbox as source text** (see ``preamble.build``): self-contained —
    ``pygenogrove`` imported inside, no module globals, no helpers from this module.
    """
    import pygenogrove as pg

    def num(v):  # a value the consensus could not give stays null, never 0
        return None if v in ("", "NA", "None", None) else int(float(v))

    created, n = [], 0
    for r in records:
        chrom, start, end = r["chrom"], int(r["start"]), int(r["end"])
        key = grove.insert(chrom, pg.GenomicCoordinate(".", start, end - 1),  # BED half-open -> closed
                           {"type": "cn_segment", "source": "PCAWG-consensus-CN",
                            "sample": r["sample"], "cohort": r.get("cohort"),
                            "total_cn": num(r["total_cn"]), "major_cn": num(r.get("major_cn")),
                            "minor_cn": num(r.get("minor_cn")), "star": num(r.get("star"))})
        created.append(("key", chrom, key))
        n += 1
    return n, created


def detach(grove, created) -> None:
    """Remove exactly the segment keys this attachment inserted."""
    for _, chrom, key in created:
        grove.remove_key(chrom, key)


LAYER = Layer(
    name="cn",
    axis="genomic",
    kind="node",
    title="Copy number — PCAWG consensus segments per tumour, per tumour cohort",
    when="a question asks whether a gene or region is amplified, deleted or copy-neutral in "
         "tumours, or needs copy number as a confounder — is an expression gain explained by a "
         "gain of the gene — in one or more PCAWG tumour cohorts; the sandbox filters on CN_COHORTS",
    schema='indexed `{"type":"cn_segment", "source":"PCAWG-consensus-CN", "sample":<WGS aliquot '
           'id>, "cohort":<PCAWG project code>, "total_cn":<int|null>, "major_cn":<int|null>, '
           '"minor_cn":<int|null>, "star":<1|2|3|null, consensus confidence>}` segment nodes '
           '(a null is a segment the consensus could not call — skip it, never treat it as 0), '
           'returned by `intersect` like genes and cCREs — filter on type, then on sample. The '
           'copy number of a gene in a tumour is the segment(s) under its span for that sample; '
           'a gene spanning two segments has two values, report both. Keep `cohort in '
           'CN_COHORTS`. No edge joins a segment to a gene: overlap is always a spatial query',
    attach=attach,
)
