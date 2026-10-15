# SPDX-License-Identifier: GPL-3.0-or-later
"""Command-line entry point for canopy.

A thin wrapper: parse the question, hand the agent loop (:mod:`genogrove_canopy.llm`) a
:class:`Session` — the host side of its tools, which attach layers into the grove, build own
data into it and run programs under restrictions (:mod:`genogrove_canopy.sandbox`) — and print
the answer. The host
resolves each dataset to a serialized ``.gg`` and binds it as ``GROVE``; the generated code
only queries.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from functools import lru_cache
from pathlib import Path

from genogrove_canopy import __version__, llm, log, preamble, resources, sandbox

# Default Anthropic model for code generation. Opus is the most capable tier and
# the connected-interval reasoning here is the paper's headline contribution, so
# we do not downgrade by default.
DEFAULT_MODEL = "claude-opus-4-8"

# The shipped grove: GENCODE structure + the ENCODE cCRE layer, in one pinned artifact.
# The enhancer→gene layer is *not* in it — it is cohort-specific and resolved per question
# (see ``layers.enhancers`` and ``_answer``).
_BASE = "gencode.human"

# When a question needs enhancers but names no tissue, load this cohort and say so.
DEFAULT_COHORT = "EFO:0005726"  # LNCaP clone FGC (prostate cancer) — the flagship cohort

#: The read-only ``GroveView`` methods advertised to the model, in the order the prompt lists them.
#:
#: A single object rather than prose so the contract cannot drift from the pinned build:
#: ``tests/test_api_surface.py`` asserts every name here exists on ``pg.GroveView``. Advertising a
#: method the build lacks is the failure that matters — generated code raises ``AttributeError``
#: inside the sandbox and the user sees a broken answer, not a build error. It has happened in both
#: directions: 457feaa removed ``get_edge_list`` as absent, and it is present on 0.7.4.
QUERY_SURFACE = ("intersect", "flanking", "get_neighbors", "get_edges", "get_edge_list",
                 "get_neighbors_if")


def build_parser() -> argparse.ArgumentParser:
    """Construct the argument parser for the ``canopy`` command."""
    parser = argparse.ArgumentParser(
        prog="canopy",
        description="Ask plain-English questions over connected genomic intervals.",
    )
    parser.add_argument("question", nargs="?", help="The natural-language question to answer.")
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Anthropic model to use for code generation (default: {DEFAULT_MODEL}).",
    )
    parser.add_argument(
        "--format",
        choices=("text", "bed", "tsv", "json"),
        default="text",
        help="Output format (default: text — an aligned table showing every field, led by "
             "the summary line; bed/tsv/json are machine formats). Scalar answers ignore this.",
    )
    parser.add_argument(
        "--show-code",
        action="store_true",
        help="Print each program the agent runs, and what came back, as it goes (stderr).",
    )
    parser.add_argument(
        "--trace",
        metavar="FILE",
        help="Write the full trace — every tool call, its input and output — as JSON to FILE.",
    )
    parser.add_argument(
        "--init",
        action="store_true",
        help="Download the shipped grove now (the pinned ~109 MB .gg) and exit, so the "
             "first real query is instant.",
    )
    parser.add_argument(
        "-i", "--interactive",
        action="store_true",
        help="Interactive session: keep the grove(s) open across questions (the ~200 ms "
             "open is paid once, then queries are sub-ms). One question per line.",
    )
    parser.add_argument(
        "--cohort",
        action="append",
        metavar="NAME",
        help="Load the enhancer→gene layer for this ENCODE-rE2G cohort (name or ontology id, "
             "e.g. 'LNCaP' or 'EFO:0005726'; repeatable for several). Enables enhancer "
             "queries. Omit and enhancers still load for an enhancer question, defaulting to "
             f"the flagship cohort. See --list-cohorts.",
    )
    parser.add_argument(
        "--data",
        action="append",
        metavar="FILE",
        default=[],
        help="Your own interval file (BED, TSV, …; repeatable). The agent may read it and put "
             "its records into the grove with `build`, then query them like any layer. Only "
             "these files are readable beyond the pinned datasets.",
    )
    parser.add_argument(
        "--list-cohorts",
        action="store_true",
        help="List the available ENCODE-rE2G cohorts (name, ontology id, type, replicates) "
             "and exit.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


_BRIDGE = Path(__file__).parent / "data" / "cohorts.tsv"


@lru_cache(maxsize=1)
def _bridge() -> list[dict]:
    """The curated tissue-term rows of ``data/cohorts.tsv`` (see its header)."""
    import csv

    with _BRIDGE.open() as fh:
        rows = list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))
    for r in rows:
        r["aliases"] = [a for a in r["aliases"].split(";") if a]
        r["re2g"] = [r["re2g"]] if r["re2g"] else []
        r["pcawg"] = [c for c in r["pcawg"].split(",") if c]
        r["hichip"] = [c for c in r["hichip"].split(",") if c]
    return rows


def _resolve_cohorts(specs):
    """Map cohort specs (``--cohort`` values or the model's ``attach`` cohort terms) to per-layer keys:
    ``{label: {"re2g": [ontology ids], "pcawg": [project codes], "hichip": [TCGA projects]}}``.

    Each spec resolves, in order, as: a layer key as written (an rE2G ontology id, a PCAWG
    project code, a TCGA project); an exact rE2G biosample name (a bridge row for the same word
    adds its PCAWG codes and TCGA projects); a bridge term or alias (``data/cohorts.tsv`` — one
    tissue word gives every layer's keys, which is what lets one cohort term drive
    enhancers, SVs *and* contacts); or a case-insensitive substring of an rE2G biosample name
    (most-replicated match wins). Raises ``SystemExit`` with a pointer to ``--list-cohorts`` on
    no match.
    """
    catalog = resources.re2g_cohorts()
    pcawg = {r["project_code"] for r in resources.pcawg_cohorts()}
    hichip = {r["project_code"] for r in resources.hichip_cohorts()}
    none = {"re2g": [], "pcawg": [], "hichip": []}
    chosen = {}
    for spec in specs:
        s = spec.strip().lower()
        if not s:
            raise SystemExit("canopy: empty cohort — see --list-cohorts")
        hit = next((c for c in catalog if c["ontology_id"].lower() == s), None)
        if hit:
            chosen[hit["name"]] = {**none, "re2g": [hit["ontology_id"]]}
            continue
        if spec.strip().upper() in pcawg:
            chosen[spec.strip().upper()] = {**none, "pcawg": [spec.strip().upper()]}
            continue
        if spec.strip().upper() in hichip:
            chosen[spec.strip().upper()] = {**none, "hichip": [spec.strip().upper()]}
            continue
        exact = next((c for c in catalog if c["name"].strip().lower() == s), None)
        row = next((r for r in _bridge() if s == r["term"] or s in (a.lower() for a in r["aliases"])), None)
        if exact:  # a bridge row for the same word names this very biosample (tested), so its
            chosen[exact["name"]] = {"re2g": [exact["ontology_id"]],  # tumour cohorts come along
                                     "pcawg": row["pcawg"] if row else [],
                                     "hichip": row["hichip"] if row else []}
            continue
        if row:
            chosen[row["term"]] = {"re2g": row["re2g"], "pcawg": row["pcawg"], "hichip": row["hichip"]}
            continue
        hit = next((c for c in catalog if s in c["name"].lower()), None)
        if hit is None:
            raise SystemExit(f"canopy: no cohort matches {spec!r} — see --list-cohorts")
        chosen[hit["name"]] = {**none, "re2g": [hit["ontology_id"]]}
    return chosen


def _list_cohorts() -> None:
    """Print what ``--cohort`` / ``attach`` can name: the tissue terms that drive both layers,
    then each layer's own keys (rE2G biosamples most-replicated first, PCAWG project codes)."""
    print("# tissue terms (data/cohorts.tsv) — one word resolves every layer")
    print(f"{'term':14}  {'rE2G':14}  {'PCAWG':32}  {'HiChIP':20}  aliases")
    for r in _bridge():
        print(f"{r['term']:14}  {(r['re2g'] or ['-'])[0]:14}  {','.join(r['pcawg']) or '-':32}  "
              f"{','.join(r['hichip']) or '-':20}  {'; '.join(r['aliases'])}")
    print("\n# ENCODE-rE2G biosamples (enhancers)")
    print(f"{'ontology id':16}  {'reps':>4}  {'type':16}  name")
    for c in resources.re2g_cohorts():
        print(f"{c['ontology_id']:16}  {c['n_replicates']:>4}  {c['type']:16}  {c['name']}")
    print("\n# PCAWG project codes (structural variants)")
    print(f"{'code':10}  {'samples':>7}  {'SVs':>7}")
    for r in resources.pcawg_cohorts():
        print(f"{r['project_code']:10}  {r['n_samples']:>7}  {r['n_svs']:>7}")
    print("\n# TCGA projects (H3K27ac HiChIP contacts)")
    print(f"{'project':10}  {'samples':>7}  {'loops':>8}")
    for r in resources.hichip_cohorts():
        print(f"{r['project_code']:10}  {r['n_samples']:>7}  {r['n_loops']:>8}")


def _grove_context(own=()):
    """Resolve the shipped grove to ``(resources_block, code_preamble, data_paths)``. ``own`` are
    the user's ``--data`` files: granted to the sandbox and named in the block for ``build``.

    The grove is the GENCODE backbone with the **Tier-1 static layers already in it** — currently
    the ENCODE cCRE registry, built into the pinned artifact rather than baked on first run (see
    ``resources.ensure_all_grove``). The preamble binds ``GROVE`` to an open ``GroveView`` of it,
    so one `intersect` returns genes *and* cCREs. The enhancer layer is **not** in the artifact
    (it is cohort-specific): once the model calls ``attach``, ``Session`` runs
    ``preamble.build(gg, cohort_links, sv_files)`` **instead of** this base preamble — the two
    are exclusive, never concatenated: the base one evicts ``_CANOPY_STATE`` from the query's
    namespace, the cohort one memoises the attached grove in it. Both bind ``GROVE``, so
    generated code never opens a path itself.
    """
    from genogrove_canopy import layers
    from genogrove_canopy.layers import enhancers, hichip, sv

    # Resolved: the sandbox compares every read against `Path.resolve()`d roots, so a symlinked
    # cache dir spelled two ways would refuse its own grove.
    gg = str(resources.ensure_all_grove(_BASE).resolve())
    # Landlock grants directories by inode; create layer roots before the worker starts.
    enhancers.LINKS_DIR.mkdir(parents=True, exist_ok=True)
    sv.SV_DIR.mkdir(parents=True, exist_ok=True)
    hichip.HICHIP_DIR.mkdir(parents=True, exist_ok=True)
    own = [str(Path(p).resolve()) for p in own]
    for p in own:
        if not Path(p).is_file():
            raise SystemExit(f"canopy: --data {p}: not a file")
    block = resources_block(
        "GROVE", resources.RESOURCES[_BASE].description,
        layers.catalogue_block(["ccre", "enhancers", "sv", "hichip"]), own,
    )
    # The sandbox reads only these roots. `LINKS_DIR` is where `enhancers.preamble`'s
    # `attach_links` opens a cohort's links table, so it must be granted alongside the grove.
    return block, preamble.build(gg), [gg, str(enhancers.LINKS_DIR.resolve()), str(sv.SV_DIR.resolve()),
                                       str(hichip.HICHIP_DIR.resolve()), *own]


def resources_block(var: str, description: str, layers_block: str, own=()) -> str:
    """Render the "Available resources" text the model reads, for grove handle ``var``.

    Pure so it can be tested without a grove: ``_grove_context`` resolves the ~109 MB artifact,
    which a unit test has no business downloading. That matters here because this is the only place
    ``QUERY_SURFACE`` reaches the model — a test asserting those names appear in the *rendered*
    block is checking the real contract, not that a token exists somewhere in a source file.
    """
    return (
        f"- `{var}`: an **open** grove handle ({description}) — gene/transcript/exon structure "
        f"**plus the ENCODE cCRE nodes**, and, when you call `attach` (see "
        f"\"Per-question layers\"), that cohort's rE2G enhancer nodes/edges, PCAWG "
        f"breakpoint edges and/or TCGA HiChIP contact windows and edges, attached by the host "
        f"before your next `run`. Query `{var}` directly. "
        f"**Never open a path yourself** — a handle you "
        f"open lacks the attached layer. A **located** query (a variant at chr7:55191822) reads "
        f"just that locus; a **genome-wide / gene-name** query works from the same handle. "
        f"Read-only in `run` — mutators raise; query with: {', '.join(f'`{m}`' for m in QUERY_SURFACE)}. "
        f"In `build` it also takes `insert`, `add_edge`, `add_external_key`, nothing else.\n"
        f"  Layers in the grove — nodes come back from `intersect` alongside genes, filter on "
        f"`source`/`type`:\n"
        f"  {layers_block.replace(chr(10), chr(10) + '  ')}\n"
        f"- `attach` cohort terms that resolve every layer (data/cohorts.tsv): "
        f"{', '.join(r['term'] for r in _bridge())}. A plain tissue word means the tissue "
        f"biosample for enhancers; the disease word ('liver cancer', 'HCC') means the cancer cell "
        f"line. An ENCODE biosample name/id, a PCAWG project code or a TCGA project (TCGA-BRCA) "
        f"resolves one layer directly.\n"
        + ("".join(f"- Your own file, readable with `open()`, to put into `{var}` with `build` (see "
                   f"\"Own data\"): `{p}`\n" for p in own) if own else
           "- No own files were given (`--data FILE` adds one), so `build` has nothing to read.\n")
    )


def _parse_output(text: str):
    """Split the generated code's stdout into ``(records, passthrough)``.

    JSONL dict lines are the feature records; every other non-empty line is the agent's
    ``label: value`` summary. Shared by the text renderer and the ``serve`` web layer (which
    wants the records as structured data to draw, not a pre-formatted table)."""
    records, passthrough = [], []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        try:
            obj = json.loads(s)
        except ValueError:
            obj = None
        if isinstance(obj, dict):
            records.append(obj)
        else:
            passthrough.append(line)  # a summary / scalar line — shown before the table
    return records, passthrough


def _render(text: str, fmt: str) -> str:
    """Render the generated code's stdout. Non-JSON lines (the agent's ``label: value``
    summary) **lead**, then the JSONL feature records become the chosen ``fmt`` table."""
    records, passthrough = _parse_output(text)
    out = list(passthrough)
    if records:
        out.append(_format_records(records, fmt))
    return "\n".join(p for p in out if p) + "\n"


def _record_columns(records: list[dict]) -> list[str]:
    """Column order across (possibly heterogeneous) records: coordinates, then identity,
    then edge evidence, then anything else — union of keys, stable order."""
    order = ["chrom", "start", "end", "strand", "name", "type", "class",
             "score", "n", "cohort", "target", "id", "biotype"]
    cols = [k for k in order if any(k in r for r in records)]
    for r in records:  # append any keys the agent used that aren't in the preferred order
        for k in r:
            if k not in cols:
                cols.append(k)
    return cols



#: Columns dropped from the text table when they carry no per-row information. Not applied to
#: tsv/json/bed — a machine format must stay a stable rectangle, and a downstream parser cannot
#: know a column vanished because it was constant.
_MIN_ROWS_TO_COLLAPSE = 2

#: Never hoisted, however constant. Every other collapsible column is metadata *about the query*
#: (`cohort`, `target`, `type`, `n`, `strand`); `chrom` is part of each row's identity. Hoisting it
#: makes a row stop being self-describing — coordinates copied out of the table without their
#: chromosome are wrong, not merely incomplete. Worth the ~7 characters.
_NEVER_COLLAPSED = frozenset({"chrom"})


def _constant_columns(records: list[dict], cols: list[str]) -> dict:
    """Columns whose value is identical in every record, as ``{column: value}``.

    An enhancer answer repeats `type=enhancer`, `cohort=EFO:0005726`, `target=AR` and a `.` strand
    down 22 rows — five columns of no information, pushing the ones that differ off the screen.
    Stated once above the table instead.
    """
    if len(records) < _MIN_ROWS_TO_COLLAPSE:
        return {}
    shared = {}
    for c in cols:
        if c in _NEVER_COLLAPSED:
            continue
        values = {_cell(r.get(c, "")) for r in records}
        if len(values) == 1 and all(c in r for r in records):
            value = values.pop()
            if value != "":
                shared[c] = value
    # Never collapse everything away: a table with no columns is not a table.
    return shared if len(shared) < len(cols) else {}


def _cell(value) -> str:
    """Render one cell for the text table.

    A list of dicts — `ccre_overlap` is the one in practice — becomes `pELS:332 PLS:277 PLS:162`
    rather than 200 characters of `repr`. The full structure stays in `--format json`, which is
    where a caller who wants the ids should be looking.
    """
    if isinstance(value, list) and value and all(isinstance(v, dict) for v in value):
        parts = []
        for item in value:
            label = item.get("class") or item.get("type") or item.get("id") or "?"
            size = item.get("bp")
            parts.append(f"{label}:{size}" if size is not None else str(label))
        return " ".join(parts)
    if isinstance(value, list) and not value:
        return "-"  # an empty list is a finding (no cCRE overlap), not a missing value
    if value is None:
        return "."  # a null field, in the BED/VCF missing-value convention (#37)
    return str(value)


def _format_records(records: list[dict], fmt: str) -> str:
    if fmt == "json":
        return "\n".join(json.dumps(r) for r in records)
    cols = _record_columns(records)
    if fmt in ("text", "tsv"):
        if fmt == "tsv":
            rows = ["\t".join(cols)]
            # null → "." (#37); an absent key stays "" (unset, not set-to-null)
            rows += ["\t".join("." if r.get(c, "") is None else str(r.get(c, "")) for c in cols)
                     for r in records]
            return "\n".join(rows)
        # text: an aligned table for a human to read. Two things keep it readable that the
        # machine formats deliberately skip — see `_constant_columns` and `_cell`.
        shared = _constant_columns(records, cols)
        cols = [c for c in cols if c not in shared]
        body = [{c: _cell(r.get(c, "")) for c in cols} for r in records]
        width = {c: max(len(c), max((len(b[c]) for b in body), default=0)) for c in cols}
        fmt_row = lambda vals: "  ".join(str(v).ljust(width[c]) for c, v in zip(cols, vals)).rstrip()
        rows = []
        if shared:  # stated once instead of repeated down every row
            rows.append("  ".join(f"{k}={v}" for k, v in shared.items()))
        rows.append(fmt_row(cols))
        rows += [fmt_row([b[c] for c in cols]) for b in body]
        return "\n".join(rows)
    # BED: 0-based closed -> half-open (end + 1); host owns the conversion, once.
    rows = ["#chrom\tstart\tend\tname\tscore\tstrand"]
    for r in records:
        if "start" not in r or "end" not in r:
            rows.append(json.dumps(r))  # not an interval record; emit verbatim
            continue
        rows.append("\t".join(str(v) for v in (
            r.get("chrom", "."), r["start"], int(r["end"]) + 1,
            r.get("name") or r.get("id") or ".",
            "." if r.get("score") is None else r["score"],    # absent or null → "." (#37)
            "." if r.get("strand") is None else r["strand"],
        )))
    return "\n".join(rows)


def _pygenogrove_site_dir() -> str:
    """The site-packages dir holding ``pygenogrove``, for the sandbox's sys.path."""
    import pygenogrove

    f = Path(pygenogrove.__file__).resolve()
    return str(f.parent.parent if f.name == "__init__.py" else f.parent)


def _resolve_query_cohorts(override, cohort_hint, layers):
    """Grounded cohort resolution for an ``attach`` of ``layers``, precedence-ordered: the
    user's ``override`` specs (``--cohort`` / the picker) → the model's ``cohort`` (resolved
    against the catalogs; no match → none, so a wrong tissue is never silently substituted) →
    the default cohort, but **only for an enhancer question**: SVs are per tumour cohort and
    there is no honest default. Returns ``(cohorts, note)`` — the ``_resolve_cohorts`` map and
    a note or ``None``."""
    if override:
        return _resolve_cohorts(override), None
    if cohort_hint:  # one or more, `;`-separated — a comparison question names several
        try:
            return _resolve_cohorts([c for c in map(str.strip, cohort_hint.split(";")) if c]), None
        except SystemExit:  # the model named a tissue with no catalog match — don't substitute
            return {}, f"no cohort matched {cohort_hint!r} — nothing attached (see --list-cohorts)"
    if "enhancers" in layers:
        return _resolve_cohorts([DEFAULT_COHORT]), "default"
    return {}, "SVs and contacts are per tumour cohort — name a tissue or pass --cohort; nothing attached"


def prepare_layers(cohorts, layers, say):
    """Materialise the declared ``layers`` for the resolved ``cohorts`` (host side, shared by
    the CLI and ``serve``): returns ``(cohort_links, sv_files, hichip_files)`` for
    ``preamble.build``. ``say(text)`` reports progress and, per declared layer that has **no**
    key in these cohorts, says so — the generated code would otherwise run with an empty
    ``COHORTS`` / ``SV_COHORTS`` / ``HICHIP_COHORTS`` and print a plausible-looking zero."""
    from genogrove_canopy.layers import enhancers, hichip, sv

    cohort_links, sv_files, hichip_files = {}, {}, {}
    names = "; ".join(cohorts)
    if "enhancers" in layers:
        ids = _cohort_ids(cohorts)
        if ids:
            say(f"Loading ENCODE-rE2G links — cohort(s) {names}")
            cohort_links = {cid: str(enhancers.links_file(cid)) for cid in ids if enhancers.ensure_index(cid)}
            if not cohort_links:
                say(f"no rE2G index for cohort(s) {names} — no enhancers attached")
        elif cohorts:
            say(f"no rE2G biosample for cohort(s) {names} — no enhancers attached")
    if "sv" in layers:
        codes = _pcawg_codes(cohorts)
        if codes:
            say(f"Loading PCAWG SVs — cohort(s) {'; '.join(codes)}")
            sv_files = {c: str(sv.cohort_file(c)) for c in codes}
        elif cohorts:
            say(f"no PCAWG cohort for cohort(s) {names} — no SVs attached")
    if "hichip" in layers:
        codes = _hichip_codes(cohorts)
        if codes:
            say(f"Loading TCGA HiChIP loops — cohort(s) {'; '.join(codes)}")
            hichip_files = {c: str(hichip.cohort_file(c)) for c in codes}
        elif cohorts:
            say(f"no TCGA HiChIP cohort for cohort(s) {names} — no contacts attached")
    return cohort_links, sv_files, hichip_files


def _prepare() -> None:
    """First-run notice if the shipped grove isn't cached yet: a one-time ~109 MB download of the
    pinned unified `.gg` (gene structure + cCREs). No local build."""
    if not resources._all_grove_gg(_BASE).exists():
        log.say(f"Fetching the grove ({resources.RESOURCES[_BASE].grove_contents}) — "
                "first run only, a pinned ~109 MB .gg")


def _cohort_ids(cohorts) -> list:
    """The rE2G ontology ids across the resolved ``cohorts`` (the enhancer index is keyed by id)."""
    return [i for c in cohorts.values() for i in c["re2g"]]


def _pcawg_codes(cohorts) -> list:
    """The PCAWG project codes across the resolved ``cohorts`` (the SV cohort unit)."""
    return [i for c in cohorts.values() for i in c["pcawg"]]


def _hichip_codes(cohorts) -> list:
    """The TCGA projects across the resolved ``cohorts`` (the HiChIP cohort unit)."""
    return [i for c in cohorts.values() for i in c["hichip"]]



class Session:
    """The host side of the agent's three tools, for one question over one grove.

    ``attach(layers, cohort)`` grounds the cohort(s) (``override`` wins — validated by the caller
    before the loop starts — then the model's term, then the enhancer default), materialises the
    layers' data and rebuilds the preamble so the next ``run`` sees them; ``build(code)`` runs a
    program with the writers enabled; ``run(code)`` executes one program through
    ``execute(script)`` (a ``Worker.submit``) and returns its stdout, or the error, as the text
    the model reads. A bad call comes back as text; anything that does raise (a download
    failing) is turned into an error result by the loop, so the model is the one that reacts.

    Attached layers accumulate across calls within the question, so ``COHORTS`` names every
    cohort the model asked for. Across questions the sandbox's warm grove keeps what was
    attached (``_CANOPY_STATE``), so attaching again is a memo hit.
    """

    def __init__(self, *, gg, base, execute, override=(), say=log.say, warm=False):
        self.gg, self.base, self.execute, self.override, self.say = gg, base, execute, list(override), say
        self.cohort_links, self.sv_files, self.hichip_files = {}, {}, {}
        self.notes = []         # what was attached / what was missing, for the UI
        self.stdout = ""        # the last successful run's output — the rows under the answer
        self.code = ""          # the last program run
        # Every run uses the warm (memoised) grove, not the lazy view: after a `build` in this
        # question, or from the start when own files are given — a build in an earlier question
        # of the session may have put them in, and the lazy view would not see that.
        self.built = warm

    def attach(self, layers, cohort=""):
        if not layers:
            return "error: name at least one layer (enhancers, sv, hichip)"
        cohorts, note = _resolve_query_cohorts(self.override, cohort, layers)
        said = []
        if note and note != "default":  # the reason nothing resolves: the model must read it too
            self.say(note)
            self.notes.append(note)
            said.append(note)
        links, svs, loops = prepare_layers(cohorts, layers, lambda m: (self.say(m), said.append(m)))
        self.notes += [m for m in said if m.startswith("no ")]
        self.cohort_links.update(links)
        self.sv_files.update(svs)
        self.hichip_files.update(loops)
        attached = []
        if links or svs or loops:
            what = " + ".join(w for w, d in (("enhancers", links), ("SVs", svs), ("contacts", loops)) if d)
            attached.append(f"{what} attached from {'; '.join(cohorts)}"
                            + (" (default — the question named no tissue)" if note == "default" else ""))
            self.notes.append(attached[-1])
        return "\n".join(attached + [m for m in said if m not in attached and not m.startswith("Loading")] + [
            f"COHORTS = {json.dumps(sorted(self.cohort_links))}",
            f"SV_COHORTS = {json.dumps(sorted(self.sv_files))}",
            f"HICHIP_COHORTS = {json.dumps(sorted(self.hichip_files))}",
        ])

    def _preamble(self, mutable=False) -> str:
        if self.cohort_links or self.sv_files or self.hichip_files or self.built or mutable:
            return preamble.build(self.gg, self.cohort_links, self.sv_files, self.hichip_files,
                                  mutable=mutable, warm=self.built)
        return self.base

    def _execute(self, code, mutable=False):
        # JSONL is the output contract, so guarantee `json` is importable even if the
        # generated code forgets the import (it's already in the allowlist).
        result = self.execute("import json\n" + self._preamble(mutable) + code)
        error = sandbox.result_error(result)
        return ("error:\n" + error) if error else result.stdout

    def build(self, code):
        """The model's own data goes into the warm grove; from here on every run uses it."""
        out = self._execute(code, mutable=True)
        self.built = True   # even on failure: partial inserts are in the memo now
        self.notes.append("own data built into the grove" if not out.startswith("error:") else "a build failed")
        return out if out.strip() else "(built; the program printed nothing)"

    def run(self, code):
        out = self._execute(code)
        if out.startswith("error:"):
            return out
        # every success replaces the rows, an empty one included — the answer must not sit on
        # the rows of an earlier, superseded run
        self.code, self.stdout = code, out
        return out if out.strip() else "(the program produced no output)"


def _show_turn(entry) -> None:
    """``--show-code``: each tool call and its result, as they happen."""
    if entry["tool"] in ("run", "build"):
        print(f"# --- {entry['tool']} ---", file=sys.stderr)
        print(entry["input"].get("code", ""), file=sys.stderr)
    else:
        print(f"# --- attach {entry['input']} ---", file=sys.stderr)
    print("# --- result ---", file=sys.stderr)
    print(entry["output"], file=sys.stderr)


def _answer(question, *, system_prompt, base, gg, args, execute, warm=False):
    """Answer one question: the agent loop over a :class:`Session`, then render.

    ``execute`` is a ``script -> SandboxResult`` callable (a ``Worker.submit``). Returns
    ``(out, error, trace, seconds, built)`` — ``out`` is what goes to stdout: the model's answer
    with the last run's rows rendered in ``args.format`` beneath it for ``text``; for the machine
    formats only the rows, the answer goes to stderr so stdout stays parseable. Exactly one of
    out/error is non-empty. ``built`` says a ``build`` happened, so a REPL can keep later
    questions warm (``warm``).
    """
    session = Session(gg=gg, base=base, execute=execute, override=args.cohort or (),
                      warm=warm or bool(getattr(args, "data", ())))
    log.say(f"Asking {args.model}")
    t0 = time.perf_counter()
    text, trace = llm.answer(question, system_prompt, session, model=args.model,
                             on_turn=_show_turn if args.show_code else None)
    took = time.perf_counter() - t0
    if args.trace:
        Path(args.trace).write_text(json.dumps({"question": question, "answer": text, "trace": trace}, indent=1))
    if not text:
        return "", (f"no answer after {len(trace)} tool calls — the turn cap ({llm.MAX_TURNS}) "
                    "was reached or the reply was cut off"), trace, took, session.built
    rows = _render(session.stdout, args.format) if session.stdout.strip() else ""
    if args.format != "text":
        print(text, file=sys.stderr)
        return rows, "", trace, took, session.built
    return text + "\n" + (f"\n{rows}" if rows.strip() else ""), "", trace, took, session.built


def _interactive(args, *, system_prompt, preamble, data_paths, site_dir) -> int:
    """Warm-worker REPL: open the grove(s) once, then answer questions until EOF/'exit'."""
    worker = sandbox.Worker(data_paths=data_paths, extra_syspath=[site_dir])
    print("canopy interactive — one question per line; Ctrl-D or 'exit' to quit.",
          file=sys.stderr)
    warm = False  # once a question builds into the grove, every later one must use the warm grove
    try:
        while True:
            try:
                question = input("ask> ").strip()
            except (EOFError, KeyboardInterrupt):
                print(file=sys.stderr)
                break
            if not question:
                continue
            if question in ("exit", "quit"):
                break
            try:
                out, err, trace, took, built = _answer(question, system_prompt=system_prompt,
                                                       base=preamble, gg=data_paths[0], args=args,
                                                       execute=worker.submit, warm=warm)
            except Exception as exc:  # e.g. an LLM error — keep the session alive
                print(f"canopy: {exc}", file=sys.stderr)
                continue
            warm = warm or built
            if err:
                print(err, file=sys.stderr)
            else:
                sys.stdout.write(out)
                sys.stdout.flush()
            log.say(f"Answered in {took:.2f}s over {len(trace)} tool call(s)")
    finally:
        worker.close()
    return 0


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and run the end-to-end loop. Returns a process exit code."""
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "serve":  # `canopy serve` — local web front-end over the same pipeline
        from genogrove_canopy import serve
        return serve.main(argv[1:])

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list_cohorts:
        _list_cohorts()
        return 0

    if args.init:  # prime the shipped grove (GENCODE + cCREs) ahead of first use, then exit
        try:
            t0 = time.perf_counter()
            _prepare()
            resources.ensure_all_grove(_BASE)
            log.took("Ready", time.perf_counter() - t0)
        except Exception as exc:
            print(f"canopy: {exc}", file=sys.stderr)
            return 1
        return 0

    if not args.question and not args.interactive:
        parser.print_help()
        return 0

    try:
        if args.cohort:  # a bad override must fail here, before any model turn is paid for
            _resolve_cohorts(args.cohort)
        # The grove is cohort-independent (GENCODE + cCREs). Layers are attached per
        # question when the model calls `attach` — see Session.
        _prepare()
        resources_block, preamble, data_paths = _grove_context(args.data)
        site_dir = _pygenogrove_site_dir()
        system_prompt = llm.build_system_prompt(resources_block)
    except SystemExit:
        raise  # a clean --cohort resolution error already carries its message
    except Exception as exc:  # surface a clean message, not a traceback
        print(f"canopy: {exc}", file=sys.stderr)
        return 1

    if args.interactive:  # warm worker: grove open paid once for the whole session
        return _interactive(args, system_prompt=system_prompt, preamble=preamble,
                            data_paths=data_paths, site_dir=site_dir)

    # One-shot is still a warm worker: the agent may run several programs for one question,
    # and the memoised grove makes every run after the first cheap.
    try:
        with sandbox.Worker(data_paths=data_paths, extra_syspath=[site_dir]) as worker:
            out, err, _trace, _took, _built = _answer(
                args.question, system_prompt=system_prompt, base=preamble, gg=data_paths[0],
                args=args, execute=worker.submit,
            )
    except Exception as exc:
        print(f"canopy: {exc}", file=sys.stderr)
        return 1
    if err:
        print(err, file=sys.stderr)
        return 1
    sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
