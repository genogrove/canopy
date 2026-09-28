# SPDX-License-Identifier: GPL-3.0-or-later
"""HiChIP contact layer — loops as edges between indexed 10 kb windows.

Source-agnostic: needs only loop calls as bin pairs (two half-open BED intervals, a contact
count, a q-value, per-bin peak flags, per sample). The TCGA H3K27ac HiChIP FitHiChIP calls
(Chang lab, Nat Genet 2025; issue #40) are the validated example source.

A loop is a link between two **superintervals**: the window is the unit of evidence and
everything finer than its width is inference this layer must not do. FitHiChIP's windows are
a fixed genome-wide 10 kb grid — every anchor is 10 kb, starts on a multiple of 10 kb, and no
two distinct windows overlap — so a window is fully identified by ``(chrom, start)``. So:

* one **indexed** key per distinct window (``type="hichip_anchor"``), exact coordinates, never
  merged into a cCRE or gene, inserted on demand the first time any cohort mentions it and
  shared by every cohort attached afterwards (the caller-owned ``windows`` cache, same role as
  the enhancer layer's ``nodes``). Indexed rather than external so that **no overlap edges are
  stored**: gene → windows is ``intersect`` on the gene span filtered on type, window → genes
  and cCREs is ``intersect`` on the window. Both directions spatial, nothing materialised —
  measured at ~17 overlap edges per window otherwise, nearly all to cCREs. Windows tile, so a
  point query gains at most one extra hit, and generated code already filters hits on type;
* one ``contact_edge`` per loop per sample between the two windows, both directions, carrying
  sample, cohort, contact count, q-value and the per-sample peak flags (per-sample evidence
  belongs on the edge, not the node).

Windows live only in sessions that ask about contacts — never in the pinned artifact. A window
in a gene desert has only its contact edge: no stand-in bin, the window already is its own
interval node. Query from a gene: ``intersect(gene span)`` → windows → ``contact_edge`` →
partner window → ``intersect`` on the partner's coordinate.

Attached per question by **tumour cohort** (TCGA project), additively onto the warm grove,
with the same ``attach_tracked``/``detach`` bookkeeping as the SV layer: ``detach`` removes the
edges it made and any window this call inserted that no edge uses any more.
"""

from __future__ import annotations

from pathlib import Path

from genogrove_canopy import resources
from genogrove_canopy.layers._base import Layer

# One record per loop per sample: BED half-open bins, as the derived cohort table writes them.
_FIELDS = ("chrom1", "start1", "end1", "chrom2", "start2", "end2",
           "cc", "q", "peak1", "peak2", "sample", "cohort")

HICHIP_DIR = resources._CACHE / "hichip_cohorts"


def cohort_file(code: str) -> Path:
    """One TCGA project's loops as a plain TSV the sandbox can read, extracted once from the
    pinned derived tarball (``tcga.hichip.loops``: one ``<project>.tsv`` per project, built by
    ``tools/build_hichip_cohorts.py``). The sandbox is granted only the roots the host names,
    so the host extracts into ``<cache>/hichip_cohorts/<code>.tsv``. Disposable cache."""
    import os
    import shutil
    import tarfile

    dest = (HICHIP_DIR / f"{code}.tsv").resolve()  # resolved: see cli._grove_context
    if dest.exists():
        return dest
    if code not in {r["project_code"] for r in resources.hichip_cohorts()}:
        raise KeyError(f"{code!r} is not a TCGA HiChIP cohort (see --list-cohorts)")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.tmp")
    with tarfile.open(resources.resolve("tcga.hichip.loops")) as t, tmp.open("wb") as out:
        shutil.copyfileobj(t.extractfile(f"{code}.tsv"), out)
    tmp.replace(dest)  # atomic: a half-written table must never look cached
    return dest


def attach(grove, records, windows=None) -> int:
    """``Layer``-contract entry point: attach loops to a **mutable** grove, return their count.
    ``windows`` is the ``(chrom, start) -> Key`` cache shared by every cohort on this grove."""
    return attach_tracked(grove, records, windows)[0]


def attach_tracked(grove, records, windows=None):
    """Same as ``attach``, but also returns what it created: ``("edge", a, b, payload)`` per
    loop and ``("key", chrom, Key)`` per window this call inserted, for ``detach``.

    **Shipped into the sandbox as source text** (see ``preamble.build``): self-contained —
    ``pygenogrove`` imported inside, no module globals, no helpers from this module.
    """
    import pygenogrove as pg

    windows = windows if windows is not None else {}
    created = []

    def window(chrom, start, end):  # end: BED half-open
        key = windows.get((chrom, start))
        if key is None:  # one already there (an earlier attach with another cache)?
            key = next((k for k in grove.intersect(pg.GenomicCoordinate("*", start, start), chrom)
                        if k.data.get("type") == "hichip_anchor" and k.value.start == start), None)
        if key is None:  # windows tile, so (chrom, start) is the identity
            key = grove.insert(chrom, pg.GenomicCoordinate(".", start, end - 1),
                               {"type": "hichip_anchor", "source": "HiChIP", "chrom": chrom})
            created.append(("key", chrom, key))
        windows[(chrom, start)] = key
        return key

    n = 0
    for r in records:
        a = window(r["chrom1"], int(r["start1"]), int(r["end1"]))
        b = window(r["chrom2"], int(r["start2"]), int(r["end2"]))
        edge = {"rel": "contact_edge", "sample": r["sample"], "cohort": r.get("cohort"),
                "cc": int(r["cc"]), "q": float(r["q"]),
                "peak1": int(r["peak1"]), "peak2": int(r["peak2"])}
        grove.add_edge(a, b, edge)
        if a is not b:  # a loop within one window: one self-edge, not two
            grove.add_edge(b, a, edge)
        created.append(("edge", a, b, edge))
        n += 1
    return n, created


def detach(grove, created, windows=None) -> None:
    """Remove this attachment's contact edges (one occurrence per tracked directed edge, so an
    identical loop from another attachment survives), then every window its edges touched that
    no edge uses any more — whichever attachment inserted it, since only this layer makes
    ``hichip_anchor`` keys — forgetting them in ``windows`` too, so a later attach re-inserts."""
    import json
    from collections import Counter

    outgoing = {}
    touched = {}
    for entry in created:
        if entry[0] == "key":
            _, chrom, key = entry
            touched[id(key)] = (chrom, key)
            continue
        _, a, b, payload = entry
        for source, target in ([(a, b)] if a is b else [(a, b), (b, a)]):
            _, counts = outgoing.setdefault(id(source), (source, Counter()))
            counts[(id(target), json.dumps(payload, sort_keys=True))] += 1
            touched[id(source)] = (source.data["chrom"], source)

    # The pinned API cannot remove a selected parallel edge: rebuild the affected adjacency
    # lists, retaining other attachments' edges in their original order.
    for source, counts in outgoing.values():
        keep = []
        for target, payload in grove.get_edge_list(source):
            match = (id(target), json.dumps(payload, sort_keys=True))
            if payload and payload.get("rel") == "contact_edge" and counts[match]:
                counts[match] -= 1
            else:
                keep.append((target, payload))
        grove.remove_edges_from(source)
        for target, payload in keep:
            grove.add_edge(source, target, payload)
    for chrom, key in touched.values():
        if not grove.get_edge_list(key) and not grove.get_in_edge_list(key):
            grove.remove_key(chrom, key)
            if windows is not None:
                windows.pop((chrom, key.value.start), None)


LAYER = Layer(
    name="hichip",
    axis="genomic",
    kind="edge",
    title="Chromatin contacts — H3K27ac HiChIP loops between 10 kb windows, per tumour cohort",
    when="a question is about 3D contact, looping, or which distal regions physically touch a "
         "gene or enhancer in tumours — one or more TCGA projects; the sandbox filters on "
         "HICHIP_COHORTS",
    schema='indexed `{"type":"hichip_anchor", "source":.., "chrom":..}` 10 kb window nodes '
           '(returned by `intersect` like genes and cCREs — filter on type); one '
           '`{"rel":"contact_edge", "sample":.., "cohort":<TCGA project>, "cc":<contact count>, '
           '"q":<FitHiChIP q>, "peak1":0|1, "peak2":0|1}` edge per loop per sample between the '
           'two windows (both directions). From a gene: `intersect(pg.GenomicCoordinate("*", '
           'gene.value.start, gene.value.end), chrom)` filtered on `type == "hichip_anchor"` → '
           'windows → `get_edge_list` for `contact_edge` → partner windows → `intersect(pg.'
           'GenomicCoordinate("*", w.value.start, w.value.end), w.data["chrom"])` for what sits '
           'at the other end (strand "*": a "." query misses stranded genes). No edge joins a '
           'window to a gene or cCRE: overlap is always a spatial query. A window is 10 kb: '
           'report every gene/cCRE it overlaps, never pick one',
    attach=attach,
)
