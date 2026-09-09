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

The graph finds junctions involving MYC by walking from its gene node. Shared gene/bin
anchors can join unrelated breakpoints, so graph cycles or clusters do not establish
chromoplexy or chromothripsis. Genes are never cut; no derivative chromosome is stored.

``svclass`` is the source caller's classification, with h2hINV/t2tINV normalized to
INV; ``source_svclass`` retains the exact original label. For lifted PCAWG calls these
describe the hg19 call. ``junction_class`` describes the current coordinates/strands:
DEL-like, DUP-like, h2hINV, t2tINV, or TRA (INS retains its explicit call). This junction
geometry can differ after liftover and does not establish copy number or functional effect.

Attached per question by **tumour cohort** (``cohort_file`` extracts one PCAWG project
code's samples into a plain table; ``attach_tracked`` inserts its edges into an
already-deserialized grove, additively — a warm session holds every cohort asked about,
and edges carry ``cohort`` and ``sample`` so the generated code filters on ``SV_COHORTS``).
``attach_tracked`` returns exactly what it created and ``detach`` removes that and nothing
else, for a caller that needs a clean working copy.
"""

from __future__ import annotations

import gzip
from pathlib import Path

from genogrove_canopy import resources
from genogrove_canopy.layers._base import Layer

# One record per SV (BEDPE-style; any source producing this column shape works).
# INS records additionally carry "length", "insertion_class", "sequence" (may be
# None if the caller couldn't resolve it — never invented).
_FIELDS = ("chrom1", "start1", "end1", "chrom2", "start2", "end2",
           "sv_id", "pe_support", "strand1", "strand2", "svclass", "svmethod")

_BIN = 1_000_000  # grid size for the intergenic anchor — see module docstring; attach_tracked repeats it


def parse_bedpe(path) -> list[dict]:
    """SV records from a BEDPE(-like) file (gzipped or plain, tab-separated, header row
    naming the columns) — the PCAWG consensus SV format, or the cohort table `cohort_file`
    writes (the same columns plus ``sample`` and ``cohort``)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return read_table(fh)


def read_table(fh) -> list[dict]:
    """Rows of an open tab-separated file whose first line names the columns.

    **Shipped into the sandbox as source text** (with ``attach_tracked``), so it must stay
    self-contained: no ``gzip`` (not on the sandbox allowlist), no module globals.
    """
    header = fh.readline().rstrip("\n").split("\t")
    return [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]


SV_DIR = resources._CACHE / "sv_cohorts"


def cohort_file(code: str) -> Path:
    """One PCAWG cohort's SVs as a plain TSV the sandbox can read, built once if needed.

    The pinned tarballs hold one gzipped BEDPE per sample; the sandbox has no ``gzip`` and is
    granted only the roots the host names, so the host extracts the cohort's samples (from the
    packaged catalog's aliquot list) into ``<cache>/sv_cohorts/<code>.tsv``: the 12 BEDPE
    columns plus ``sample`` (aliquot id) and ``cohort`` (project code). Disposable cache.
    """
    import os
    import tarfile

    dest = (SV_DIR / f"{code}.tsv").resolve()  # resolved: see cli._grove_context
    if dest.exists():
        return dest
    row = next((r for r in resources.pcawg_cohorts() if r["project_code"] == code), None)
    if row is None:
        raise KeyError(f"{code!r} is not a PCAWG SV cohort (see --list-cohorts)")
    wanted = set(row["aliquot_ids"])
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(f"{dest.name}.{os.getpid()}.tmp")
    n = 0
    with tmp.open("w") as out:
        out.write("\t".join(_FIELDS + ("sample", "cohort")) + "\n")
        for name in ("pcawg.sv.icgc", "pcawg.sv.tcga"):
            with tarfile.open(resources.resolve(name)) as t:
                for m in t.getmembers():
                    aliquot = m.name.split("/")[-1].split(".")[0]
                    if not m.name.endswith(".bedpe.gz") or aliquot not in wanted:
                        continue
                    with gzip.open(t.extractfile(m), "rt") as fh:
                        for r in read_table(fh):
                            out.write("\t".join([r[f] for f in _FIELDS] + [aliquot, code]) + "\n")
                    n += 1
    if n != len(wanted):
        tmp.unlink()
        raise RuntimeError(f"{code}: found {n} of {len(wanted)} samples in the pinned tarballs")
    tmp.replace(dest)
    return dest


def attach(grove, records) -> int:
    """Attach one sample's SVs to a **mutable**, already-deserialized ``grove`` —
    the ``Layer``-contract entry point (``(grove, records) -> int``, matching
    ``ccres``/``enhancers``). Each record must carry a ``"sample"`` key alongside
    the ``_FIELDS`` above (and ``"cohort"``, carried onto the edge when present); one
    call may hold many samples — a whole cohort. Returns the number of SVs attached.

    This entry point doesn't track what it created — use ``attach_tracked`` when a
    caller may need ``detach`` later.

    **Shipped into the sandbox as source text** (same self-contained-function
    convention as ``enhancers.attach_links``): everything it needs is imported
    inside it or passed in, no module globals, so ``tests/`` can exercise the exact
    code the sandbox runs.
    """
    return attach_tracked(grove, records)[0]


def attach_tracked(grove, records, bin_size=1_000_000):  # literal: the def runs in the sandbox
    """Same as ``attach``, but also returns every edge and bin it created — as
    ``("edge", a, b, payload)`` and ``("key", index, Key)`` entries — so ``detach`` can remove
    exactly what this call added and nothing shared with a later sample.

    **Shipped into the sandbox as source text** (see ``preamble.build``): self-contained by
    construction — ``pygenogrove`` imported inside, the bin size a default argument rather
    than the module global, no helpers from this module.
    """
    import pygenogrove as pg

    created = []
    bins: dict[tuple, object] = {}  # (chrom, bin_start) -> bin Key, this call's own

    def anchors(chrom, pos):
        at = pg.GenomicCoordinate("*", pos, pos)
        genes = [k for k in grove.intersect(at, chrom) if k.data.get("type") == "gene"]
        if genes:
            return genes
        bin_start = (pos // bin_size) * bin_size
        key = bins.get((chrom, bin_start))
        if key is None:  # one already there (an earlier sample's, or from another bin lookup)?
            key = next((k for k in grove.intersect(at, chrom)
                        if k.data.get("type") == "intergenic_region"
                        and k.value.start == bin_start), None)
        if key is None:
            key = grove.insert(chrom, pg.GenomicCoordinate(".", bin_start, bin_start + bin_size - 1),
                               {"type": "intergenic_region", "_sv_owned": True})
            created.append(("key", chrom, key))
        bins[(chrom, bin_start)] = key
        return [key]

    n = 0
    for r in records:
        source_class = r["svclass"]
        svclass = "INV" if source_class in ("h2hINV", "t2tINV") else source_class
        if svclass == "INS":
            junction_class = "INS"  # an insertion needs caller evidence, not strand inference
        elif r["chrom1"] != r["chrom2"]:
            junction_class = "TRA"
        else:
            p1, p2 = int(r["start1"]), int(r["start2"])
            strands = r["strand1"] + r["strand2"]
            if p1 > p2:
                strands = strands[::-1]
            junction_class = {"++": "h2hINV", "--": "t2tINV"}.get(strands)
            if p1 != p2 and junction_class is None:
                junction_class = {"+-": "DEL-like", "-+": "DUP-like"}.get(strands)
        edge = {
            "rel": "breakpoint_edge", "svclass": svclass, "sv_id": r["sv_id"],
            "source_svclass": source_class, "junction_class": junction_class,
            "sample": r["sample"], "cohort": r.get("cohort"),
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
                created.append(("edge", a, b, edge))
        n += 1
    return n, created


def detach(grove, created) -> None:
    """Remove this attachment's edge multiplicities, then its unused SV-owned bins.

    The bindings decode payloads by value. Match the full payload and target, consuming one
    occurrence per tracked directed edge; identical calls retain the other copy. Bins shared
    by later attachments survive until their last edge is detached, in either detach order.
    """
    import json
    from collections import Counter

    outgoing = {}
    bins = {}
    for entry in created:
        if entry[0] == "key":
            _, chrom, key = entry
            bins[id(key)] = (chrom, key)
            continue
        _, a, b, payload = entry
        directions = [(a, b)] if a is b else [(a, b), (b, a)]
        for source, target in directions:
            _, counts = outgoing.setdefault(id(source), (source, Counter()))
            counts[(id(target), json.dumps(payload, sort_keys=True))] += 1
        for key, chrom in ((a, payload["chrom1"]), (b, payload["chrom2"])):
            if key.data.get("_sv_owned"):
                bins[id(key)] = (chrom, key)

    # The pinned API cannot remove a selected parallel edge from one source. Rebuild only
    # affected adjacency lists, retaining other attachments and their original order.
    for source, counts in outgoing.values():
        keep = []
        for target, payload in grove.get_edge_list(source):
            match = (id(target), json.dumps(payload, sort_keys=True))
            if counts[match]:
                counts[match] -= 1
            else:
                keep.append((target, payload))
        grove.remove_edges_from(source)
        for target, payload in keep:
            grove.add_edge(source, target, payload)
    for chrom, key in bins.values():
        if not grove.get_edge_list(key) and not grove.get_in_edge_list(key):
            grove.remove_key(chrom, key)


LAYER = Layer(
    name="sv",
    axis="genomic",
    kind="edge",
    title="Structural variants — PCAWG breakpoint edges, per tumour cohort",
    when="a question is about structural rearrangements, what a gene is joined to or how often "
         "it is broken across tumours, or whether a rearrangement changed a gene's neighbourhood "
         "— one or more PCAWG tumour cohorts (project codes such as BRCA-US); the sandbox "
         "filters on SV_COHORTS",
    schema='one `{"rel":"breakpoint_edge", "svclass":<DEL|DUP|INV|TRA|INS>, "sv_id":.., '
           '"source_svclass":.., "junction_class":.., "sample":<tumour aliquot id>, '
           '"cohort":<PCAWG project code>, "chrom1":.., "pos1":.., "strand1":.., "chrom2":.., "pos2":.., '
           '"strand2":.., "pe_support":..}` edge per SV (both directions) between the two '
           'backbone nodes its breakends fall in: the containing gene, or a 1 Mb '
           '`{"type":"intergenic_region"}` bin when no gene contains the position. Walk it from '
           'a gene to find junctions involving it; the exact breakpoint positions are on the '
           'edge. `svclass` is the source call with h2hINV/t2tINV normalized to INV; '
           '`source_svclass` preserves the original label (hg19 for PCAWG). '
           '`junction_class` describes current-assembly geometry: DEL-like, DUP-like, '
           'h2hINV, t2tINV, TRA, INS, or null when unresolved. Do not infer copy number '
           'or complex rearrangement diagnoses from this geometry. INS carries '
           '`"length"`/`"insertion_class"`/`"sequence"` too; the pinned PCAWG calls contain no INS',
    attach=attach,
)
