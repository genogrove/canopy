# SPDX-License-Identifier: GPL-3.0-or-later
"""ENCODE Registry of cCREs (V4) — the epigenomic **node** layer on the genomic axis.

One pinned, pre-indexed BED (ENCODE Registry of cCREs V4, GRCh38, 2,348,854 elements; Nature
2026, doi:10.1038/s41586-025-09909-9). Like the GENCODE grove it ships as a **hosted
bgzip+tabix pair** — no on-the-fly indexing on a user's machine — so ``resources.indexed_path``
downloads the ``.bed.gz`` + ``.tbi`` (sha-verified) and a query reads only its locus.

A cCRE is a pure **node**: an interval plus its class (``PLS``/``pELS``/``dELS``/``CA``/
``CA-CTCF``/``CA-H3K4me3``/``TF``/``CA-TF`` in V4). No edges — this layer broadens the genomic
axis with an epigenomic overlap target ("does this variant fall in a cCRE, and what class"),
it does not connect nodes.

``type`` is the Sequence Ontology term ``regulatory_region`` (SO:0005836) — the same vocabulary
as the backbone's ``gene``/``transcript``/``exon``, and a legal GFF3 column 3. SO has no term
for a *candidate* element, and the per-class terms all overstate: ``enhancer`` (SO:0000165) for
a dELS asserts validated function ENCODE never claims (hence "enhancer-**like** signature"),
while CA / CA-H3K4me3 / CA-CTCF aren't cis-regulatory in SO at all (``open_chromatin_region``
and ``CTCF_binding_site`` descend from ``biological_region``, not ``regulatory_region``). So the
evidence stays where GFF3 can express uncertainty — the attributes: ``class`` verbatim, and
``source == "ENCODE-SCREEN"`` marking the whole layer as candidate calls to read with care.
That also makes "any cCRE" one ``source`` test, independent of the class vocabulary.

The sandbox allowlist blocks file I/O, so the tabix read is **host-side** (like the enhancer
layer); only the resulting nodes reach the grove the generated code queries.
"""

from __future__ import annotations

import subprocess

from genogrove_canopy import resources
from genogrove_canopy.layers._base import Layer

_NAME = "encode.ccre.v4"  # RESOURCES catalog key (pinned hosted .bed.gz + .tbi)

# Columns of the ENCODE Registry cCRE BED (6-col: coordinate, the two accessions, class).
_FIELDS = ("chrom", "start", "end", "rdhs", "accession", "ccre_class")

_TYPE = "regulatory_region"  # SO:0005836 — see the module docstring on why not the per-class terms
_SOURCE = "ENCODE-SCREEN"    # GFF3 column 2; the "these are candidates" marker, class-independent


def in_region(chrom: str, start: int, end: int) -> list[dict]:
    """cCREs overlapping ``chrom:start-end`` (0-based closed — the grove convention).

    Reads only the locus through the hosted tabix index (downloaded + sha-verified on first
    use). Returns one dict per element carrying the ``_FIELDS`` columns.
    """
    path = resources.indexed_path(_NAME)  # hosted .bed.gz + .tbi pair; no local build
    region = f"{chrom}:{start + 1}-{end + 1}"  # 0-based closed -> tabix 1-based inclusive
    out = subprocess.run(["tabix", str(path), region],  # noqa: S603 — pinned local index
                         capture_output=True, text=True, check=True).stdout
    return [dict(zip(_FIELDS, ln.split("\t"))) for ln in out.splitlines() if ln]


def all_records():
    """Yield every cCRE in the registry (whole-BED stream) — the build-time input for baking the
    layer into the shipped grove (``resources.ensure_baked_grove``). Host-side only: reads the
    pinned ``.bed.gz`` directly, no tabix. ~2.35M records."""
    import gzip

    path = resources.indexed_path(_NAME)  # the hosted .bed.gz (its .tbi is unused here)
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln:
                yield dict(zip(_FIELDS, ln.split("\t")))


def attach(grove, records) -> int:
    """Insert cCRE ``records`` into a **mutable** grove as spatial nodes (intersect-queryable).

    Node layer — no edges. BED is half-open and the grove key 0-based closed, so ``end`` shifts
    by one. ``class`` is stored **verbatim** so the loader is release-agnostic across the V3/V4
    class-vocabulary change. Returns the number inserted.
    """
    # ponytail: one SO type for all classes. The per-class map (PLS->promoter, dELS->enhancer,
    # CA->open_chromatin_region, …) is a fixed function of `class` — derive it downstream if
    # something ever needs it, rather than baking an over-specific claim into 2.35M nodes.
    import pygenogrove as pg

    n = 0
    for r in records:
        grove.insert(
            r["chrom"],
            pg.GenomicCoordinate(".", int(r["start"]), int(r["end"]) - 1),
            {"type": _TYPE, "source": _SOURCE, "class": r["ccre_class"],
             "id": r["accession"], "rdhs": r["rdhs"]},
        )
        n += 1
    return n


LAYER = Layer(
    name="ccre",
    axis="genomic",
    kind="node",
    title="ENCODE Registry of cCREs (V4)",
    when="a question asks whether a locus/variant falls in a candidate cis-regulatory element "
         "(cCRE) or asks for the regulatory-element class at a position",
    schema='node `{"type":"regulatory_region", "source":"ENCODE-SCREEN", '
           '"class":<PLS|pELS|dELS|CA|CA-CTCF|CA-H3K4me3|TF|CA-TF>, '
           '"id":"EH38E…", "rdhs":"EH38D…"}` — filter cCREs on `source`, not `type`; no edges',
    attach=attach,
)
