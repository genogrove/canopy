# SPDX-License-Identifier: GPL-3.0-or-later
"""Structural-variant (SV) layer — breakpoint edges between backbone nodes.

Source-agnostic: needs only breakpoint pairs (two coordinates + strand per SV, per
sample). PCAWG consensus SV BEDPE is the validated example source, not the only one
this works with — anything producing the same column shape (``_FIELDS`` below)
attaches unmodified.

A sample's SV is **one edge between the two places its breakends fall**, nothing else
is created:

* a breakend inside a gene anchors to that gene node (every gene containing the
  position, so a breakend in an overlapping gene pair — EGFR / EGFR-AS1 — reaches both);
* a breakend outside every gene anchors to a 1 Mb ``intergenic_region`` bin for that
  position (key = ``pos // _BIN``, plain arithmetic, never derived from gene positions,
  never a nearby-but-uninvolved gene), created on demand and reused if present;
* one ``breakpoint_edge`` per SV joins the two anchors, both directions, carrying the
  exact positions and strands so nothing is lost by anchoring to a whole gene or bin.
  An insertion (``svclass="INS"``) is the same edge carrying ``length`` /
  ``insertion_class`` / ``sequence``.

So "what is now next to MYC in this sample" is one hop from the MYC node over
``breakpoint_edge``; complex rearrangements (chromoplexy, chromothripsis) are just more
edges — a chromoplexy loop is a cycle, a chromothriptic cluster is many edges piled into
one region. Genes are never cut, and no derivative-chromosome structure is stored.

Ephemeral, per sample: ``attach_tracked`` inserts into an already-deserialized grove and
returns exactly what it added (edges and any new bins); ``detach`` removes that and only
that, so one warm backbone (``sandbox.py``'s ``Worker``/``_CANOPY_STATE``) is reused
across samples without ever holding two samples' rearrangements at once.
"""

from __future__ import annotations

import gzip

from genogrove_canopy.layers._base import Layer

# One record per SV (BEDPE-style; any source producing this column shape works).
# INS records additionally carry "length", "insertion_class", "sequence" (may be
# None if the caller couldn't resolve it — never invented).
_FIELDS = ("chrom1", "start1", "end1", "chrom2", "start2", "end2",
           "sv_id", "pe_support", "strand1", "strand2", "svclass", "svmethod")

_BIN = 1_000_000  # grid size for the intergenic anchor — see module docstring


def parse_bedpe(path) -> list[dict]:
    """One sample's SV records from a BEDPE(-like) file (gzipped or plain, tab-
    separated, header row present) — matches the PCAWG consensus SV format."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        fh.readline()  # header
        return [dict(zip(_FIELDS, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]


def attach(grove, records) -> int:
    """Attach one sample's SVs to a **mutable**, already-deserialized ``grove`` —
    the ``Layer``-contract entry point (``(grove, records) -> int``, matching
    ``ccres``/``enhancers``). Each record must carry a ``"sample"`` key alongside
    the ``_FIELDS`` above; all records passed in one call are one sample's SVs.
    Returns the number of SVs attached.

    Ephemeral by construction, but this entry point doesn't track what it created
    — use ``attach_tracked`` when the grove is a warm copy that will be reused for
    a *different* sample afterward, so ``detach`` can clean up first.

    **Shipped into the sandbox as source text** (same self-contained-function
    convention as ``enhancers.attach_links``): everything it needs is imported
    inside it or passed in, no module globals, so ``tests/`` can exercise the exact
    code the sandbox runs.
    """
    return attach_tracked(grove, records)[0]


def attach_tracked(grove, records):
    """Same as ``attach``, but also returns every edge and bin it created — as
    ``("edge", a, b)`` and ``("key", index, Key)`` entries — so ``detach`` can remove
    exactly what this call added and nothing shared with a later sample.
    """
    import pygenogrove as pg

    created = []
    bins: dict[tuple, object] = {}  # (chrom, bin_start) -> bin Key, this call's own
    sample = records[0]["sample"] if records else None

    def anchors(chrom, pos):
        at = pg.GenomicCoordinate("*", pos, pos)
        genes = [k for k in grove.intersect(at, chrom) if k.data.get("type") == "gene"]
        if genes:
            return genes
        bin_start = (pos // _BIN) * _BIN
        key = bins.get((chrom, bin_start))
        if key is None:  # one already there (an earlier sample's, or from another bin lookup)?
            key = next((k for k in grove.intersect(at, chrom)
                        if k.data.get("type") == "intergenic_region"
                        and k.value.start == bin_start), None)
        if key is None:
            key = grove.insert(chrom, pg.GenomicCoordinate(".", bin_start, bin_start + _BIN - 1),
                               {"type": "intergenic_region"})
            created.append(("key", chrom, key))
        bins[(chrom, bin_start)] = key
        return [key]

    n = 0
    for r in records:
        edge = {
            "rel": "breakpoint_edge", "svclass": r["svclass"], "sv_id": r["sv_id"],
            "sample": sample,
            "chrom1": r["chrom1"], "pos1": int(r["start1"]), "strand1": r["strand1"],
            "chrom2": r["chrom2"], "pos2": int(r["start2"]), "strand2": r["strand2"],
            "pe_support": int(r["pe_support"]), "svmethod": r["svmethod"],
        }
        if r["svclass"] == "INS":
            edge["length"] = r.get("length")
            edge["insertion_class"] = r.get("insertion_class")
            edge["sequence"] = r.get("sequence")  # None if the caller couldn't resolve it
        seen = set()  # unordered anchor pairs this SV already joined: both breakends inside
        for a in anchors(r["chrom1"], int(r["start1"])):  # the same overlapping genes would
            for b in anchors(r["chrom2"], int(r["start2"])):  # otherwise yield (A,B) and (B,A)
                if frozenset((id(a), id(b))) in seen:
                    continue
                seen.add(frozenset((id(a), id(b))))
                grove.add_edge(a, b, edge)
                if a is not b:  # both breakends in one gene/bin: one self-edge, not two
                    grove.add_edge(b, a, edge)
                created.append(("edge", a, b))
        n += 1
    return n, created


def detach(grove, created) -> None:
    """Remove exactly what one ``attach_tracked`` call added: its breakpoint edges (both
    directions), then any bin nodes it created. Gene nodes are never touched.

    ponytail: ``remove_edge(a, b)`` drops the *first* edge between the pair, so this relies on
    no other layer putting an edge between two anchors — true today (backbone edges are
    gene→transcript→exon, enhancer edges are enhancer↔gene, bins have only ours). If a
    gene↔gene layer ever lands, remove by predicate on ``sample`` instead.
    """
    for entry in created:
        if entry[0] == "edge":
            _, a, b = entry
            grove.remove_edge(a, b)
            if a is not b:
                grove.remove_edge(b, a)
    for entry in created:
        if entry[0] == "key":
            _, index, key = entry
            grove.remove_key(index, key)


LAYER = Layer(
    name="sv",
    axis="genomic",
    kind="edge",
    title="Structural variants — breakpoint graph",
    when="a question is about structural rearrangement, a specific patient/sample's "
         "genome structure, chromothripsis/chromoplexy, or whether a gene's regulatory "
         "context changed due to a rearrangement — scoped to one sample at a time",
    schema='one `{"rel":"breakpoint_edge", "svclass":<DEL|DUP|INV|TRA|INS>, "sv_id":.., '
           '"sample":.., "chrom1":.., "pos1":.., "strand1":.., "chrom2":.., "pos2":.., '
           '"strand2":.., "pe_support":..}` edge per SV (both directions) between the two '
           'backbone nodes its breakends fall in: the containing gene, or a 1 Mb '
           '`{"type":"intergenic_region"}` bin when no gene contains the position. Walk it from '
           'a gene to see what this sample\'s rearrangement now puts next to it; the exact '
           'breakpoint positions are on the edge. INS carries `"length"`/`"insertion_class"`/'
           '`"sequence"` too',
    attach=attach,
)
