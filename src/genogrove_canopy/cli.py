# SPDX-License-Identifier: GPL-3.0-or-later
"""Command-line entry point for canopy.

A thin wrapper: parse the question, then orchestrate the three stages — generate
Python (:mod:`genogrove_canopy.llm`), execute it under restrictions (:mod:`genogrove_canopy.sandbox`), and
print the result. The host resolves each dataset to a serialized ``.gg`` and
injects its path as a variable; the generated code only deserializes and queries.
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
        help="Print the generated Python before running it.",
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
    return rows


def _resolve_cohorts(specs):
    """Map cohort specs (``--cohort`` values or the model's ``COHORT:`` terms) to per-layer keys:
    ``{label: {"re2g": [ontology ids], "pcawg": [project codes]}}``.

    Each spec resolves, in order, as: a layer key as written (an rE2G ontology id, a PCAWG
    project code); an exact rE2G biosample name (a bridge row for the same word adds its PCAWG
    codes); a bridge term or alias (``data/cohorts.tsv`` — one tissue word gives both layers'
    keys, which is what lets one ``COHORT:`` line drive enhancers *and* SVs); or a
    case-insensitive substring of an rE2G biosample name (most-replicated match wins). Raises
    ``SystemExit`` with a pointer to ``--list-cohorts`` on no match.
    """
    catalog = resources.re2g_cohorts()
    pcawg = {r["project_code"] for r in resources.pcawg_cohorts()}
    chosen = {}
    for spec in specs:
        s = spec.strip().lower()
        if not s:
            raise SystemExit("canopy: empty cohort — see --list-cohorts")
        hit = next((c for c in catalog if c["ontology_id"].lower() == s), None)
        if hit:
            chosen[hit["name"]] = {"re2g": [hit["ontology_id"]], "pcawg": []}
            continue
        if spec.strip().upper() in pcawg:
            chosen[spec.strip().upper()] = {"re2g": [], "pcawg": [spec.strip().upper()]}
            continue
        exact = next((c for c in catalog if c["name"].strip().lower() == s), None)
        row = next((r for r in _bridge() if s == r["term"] or s in (a.lower() for a in r["aliases"])), None)
        if exact:  # a bridge row for the same word names this very biosample (tested), so its
            chosen[exact["name"]] = {"re2g": [exact["ontology_id"]],  # PCAWG codes come along
                                     "pcawg": row["pcawg"] if row else []}
            continue
        if row:
            chosen[row["term"]] = {"re2g": row["re2g"], "pcawg": row["pcawg"]}
            continue
        hit = next((c for c in catalog if s in c["name"].lower()), None)
        if hit is None:
            raise SystemExit(f"canopy: no cohort matches {spec!r} — see --list-cohorts")
        chosen[hit["name"]] = {"re2g": [hit["ontology_id"]], "pcawg": []}
    return chosen


def _list_cohorts() -> None:
    """Print what ``--cohort`` / ``COHORT:`` can name: the tissue terms that drive both layers,
    then each layer's own keys (rE2G biosamples most-replicated first, PCAWG project codes)."""
    print("# tissue terms (data/cohorts.tsv) — one word resolves both layers")
    print(f"{'term':14}  {'rE2G':14}  {'PCAWG':32}  aliases")
    for r in _bridge():
        print(f"{r['term']:14}  {(r['re2g'] or ['-'])[0]:14}  {','.join(r['pcawg']) or '-':32}  "
              f"{'; '.join(r['aliases'])}")
    print("\n# ENCODE-rE2G biosamples (enhancers)")
    print(f"{'ontology id':16}  {'reps':>4}  {'type':16}  name")
    for c in resources.re2g_cohorts():
        print(f"{c['ontology_id']:16}  {c['n_replicates']:>4}  {c['type']:16}  {c['name']}")
    print("\n# PCAWG project codes (structural variants)")
    print(f"{'code':10}  {'samples':>7}  {'SVs':>7}")
    for r in resources.pcawg_cohorts():
        print(f"{r['project_code']:10}  {r['n_samples']:>7}  {r['n_svs']:>7}")


def _grove_context():
    """Resolve the shipped grove to ``(resources_block, code_preamble, data_paths)``.

    The grove is the GENCODE backbone with the **Tier-1 static layers already in it** — currently
    the ENCODE cCRE registry, built into the pinned artifact rather than baked on first run (see
    ``resources.ensure_all_grove``). The preamble binds ``GROVE`` to an open ``GroveView`` of it,
    so one `intersect` returns genes *and* cCREs. The enhancer layer is **not** in the artifact
    (it is cohort-specific): when the model declares ``COHORT``/``LAYERS``, ``_answer`` appends
    ``preamble.build(gg, cohort_links)``, which rebinds ``GROVE`` to a mutable copy with that
    cohort's nodes and edges attached — same name, so generated code never opens a path itself.
    """
    from genogrove_canopy import layers
    from genogrove_canopy.layers import enhancers, sv

    # Resolved: the sandbox compares every read against `Path.resolve()`d roots, so a symlinked
    # cache dir spelled two ways would refuse its own grove.
    gg = str(resources.ensure_all_grove(_BASE).resolve())
    block = resources_block(
        "GROVE", resources.RESOURCES[_BASE].description,
        layers.catalogue_block(["ccre", "enhancers", "sv"]),
    )
    # The sandbox reads only these roots. `LINKS_DIR` is where `enhancers.preamble`'s
    # `attach_links` opens a cohort's links table, so it must be granted alongside the grove.
    return block, preamble.build(gg), [gg, str(enhancers.LINKS_DIR.resolve()), str(sv.SV_DIR.resolve())]


def resources_block(var: str, description: str, layers_block: str) -> str:
    """Render the "Available resources" text the model reads, for grove handle ``var``.

    Pure so it can be tested without a grove: ``_grove_context`` resolves the ~109 MB artifact,
    which a unit test has no business downloading. That matters here because this is the only place
    ``QUERY_SURFACE`` reaches the model — a test asserting those names appear in the *rendered*
    block is checking the real contract, not that a token exists somewhere in a source file.
    """
    return (
        f"- `{var}`: an **open** grove handle ({description}) — gene/transcript/exon structure "
        f"**plus the ENCODE cCRE nodes**, and, when you declare `COHORT`/`LAYERS` (see "
        f"\"Per-question layers\"), that cohort's rE2G enhancer nodes/edges and/or PCAWG "
        f"breakpoint edges, attached by the host before your code runs. Query `{var}` directly. "
        f"**Never open a path yourself** — a handle you "
        f"open lacks the attached layer. A **located** query (a variant at chr7:55191822) reads "
        f"just that locus; a **genome-wide / gene-name** query works from the same handle. "
        f"Read-only — mutators raise; query with: {', '.join(f'`{m}`' for m in QUERY_SURFACE)}.\n"
        f"  Layers in the grove — nodes come back from `intersect` alongside genes, filter on "
        f"`source`/`type`:\n"
        f"  {layers_block.replace(chr(10), chr(10) + '  ')}\n"
        f"- `COHORT:` terms that resolve both layers (data/cohorts.tsv): "
        f"{', '.join(r['term'] for r in _bridge())}. A plain tissue word means the tissue "
        f"biosample for enhancers; the disease word ('liver cancer', 'HCC') means the cancer cell "
        f"line. An ENCODE biosample name/id or a PCAWG project code resolves one layer directly.\n"
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
    return str(value)


def _format_records(records: list[dict], fmt: str) -> str:
    if fmt == "json":
        return "\n".join(json.dumps(r) for r in records)
    cols = _record_columns(records)
    if fmt in ("text", "tsv"):
        if fmt == "tsv":
            rows = ["\t".join(cols)]
            rows += ["\t".join(str(r.get(c, "")) for c in cols) for r in records]
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
            r.get("name") or r.get("id") or ".", r.get("score", "."),
            r.get("strand", "."),
        )))
    return "\n".join(rows)


def _pygenogrove_site_dir() -> str:
    """The site-packages dir holding ``pygenogrove``, for the sandbox's sys.path."""
    import pygenogrove

    f = Path(pygenogrove.__file__).resolve()
    return str(f.parent.parent if f.name == "__init__.py" else f.parent)


def _resolve_query_cohorts(args, cohort_hint, layers):
    """Grounded cohort resolution for a query that declared ``layers``, precedence-ordered:
    ``--cohort`` override → the model's declared ``COHORT`` (resolved against the catalogs; no
    match → none, so a wrong tissue is never silently substituted) → the default cohort, but
    **only for an enhancer question**: SVs are per tumour cohort and there is no honest default.
    Returns ``(cohorts, note)`` — the ``_resolve_cohorts`` map and a stderr line or ``None``."""
    if args.cohort:
        return _resolve_cohorts(args.cohort), None
    if cohort_hint:  # one or more, `;`-separated — a comparison question names several
        try:
            return _resolve_cohorts([c for c in map(str.strip, cohort_hint.split(";")) if c]), None
        except SystemExit:  # the model named a tissue with no catalog match — don't substitute
            return {}, f"no cohort matched {cohort_hint!r} — nothing attached (see --list-cohorts)"
    if "enhancers" in layers:
        return _resolve_cohorts([DEFAULT_COHORT]), "default"
    return {}, "SVs are per tumour cohort — name a tissue or pass --cohort; nothing attached"


def prepare_layers(cohorts, layers, say):
    """Materialise the declared ``layers`` for the resolved ``cohorts`` (host side, shared by
    the CLI and ``serve``): returns ``(cohort_links, sv_files)`` for ``preamble.build``.
    ``say(text)`` reports progress and, per declared layer that has **no** key in these
    cohorts, says so — the generated code would otherwise run with an empty ``COHORTS`` /
    ``SV_COHORTS`` and print a plausible-looking zero."""
    from genogrove_canopy.layers import enhancers, sv

    cohort_links, sv_files = {}, {}
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
    return cohort_links, sv_files


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



def _answer(question, *, system_prompt, base, gg, args, execute):
    """Translate one question to code, run it via ``execute(script)``, and render.

    ``execute`` is a ``script -> SandboxResult`` callable (``sandbox.run`` for one-shot,
    ``Worker.submit`` for interactive). Returns ``(rendered_stdout, error_msg, gen_s, enh_s,
    exec_s)`` — exactly one of stdout/error is non-empty. Three times, not two: the rE2G
    attach sits between code-gen and execution, and leaving it out made the reported total
    wrong by however long it took.

    Per-question layers are resolved **per question**: the model declares ``COHORT``/``LAYERS``,
    the host grounds the cohort(s) (``--cohort`` overrides, repeatable), and each cohort's links
    are attached onto the mutable grove in the sandbox — reused warm across turns via
    ``_CANOPY_STATE`` — rather than fetched per target and injected as a list.
    """
    log.say(f"Generating a pygenogrove query ({args.model})")
    t0 = time.perf_counter()
    cohort_hint, layers, code = llm.generate_query(question, system_prompt, model=args.model)
    gen_s = time.perf_counter() - t0
    log.took("Query generated", gen_s)
    if args.show_code:
        print("# --- generated code ---", file=sys.stderr)
        print(code, file=sys.stderr)
    enh_pre, enh_s = "", 0.0
    if layers:  # a per-question layer was declared — resolve the cohort(s), prepare its data
        cohorts, note = _resolve_query_cohorts(args, cohort_hint, layers)
        if note and note != "default":
            log.say(note)
        t_enh = time.perf_counter()
        cohort_links, sv_files = prepare_layers(cohorts, layers, log.say)
        enh_s = time.perf_counter() - t_enh
        if cohort_links or sv_files:
            enh_pre = preamble.build(gg, cohort_links, sv_files)
            src = " (default — name a tissue or pass --cohort)" if note == "default" else ""
            what = " + ".join(w for w, d in (("rE2G", cohort_links), ("SV", sv_files)) if d)
            log.took(f"{what}: attached {'; '.join(cohorts)}{src}", enh_s)
    # JSONL is the output contract, so guarantee `json` is importable even if the
    # generated code forgets the import (it's already in the allowlist).
    log.say("Running the query over the grove")
    t1 = time.perf_counter()
    result = execute("import json\n" + base + enh_pre + code)
    exec_s = time.perf_counter() - t1
    error = sandbox.result_error(result)
    if error:
        return "", error, gen_s, enh_s, exec_s
    rendered = _render(result.stdout, args.format)
    if not rendered.strip():
        return "", "(the generated code produced no output)", gen_s, enh_s, exec_s
    return rendered, "", gen_s, enh_s, exec_s


def _interactive(args, *, system_prompt, preamble, data_paths, site_dir) -> int:
    """Warm-worker REPL: open the grove(s) once, then answer questions until EOF/'exit'."""
    worker = sandbox.Worker(data_paths=data_paths, extra_syspath=[site_dir])
    print("canopy interactive — one question per line; Ctrl-D or 'exit' to quit.",
          file=sys.stderr)
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
                out, err, gen_s, enh_s, exec_s = _answer(question, system_prompt=system_prompt,
                                                  base=preamble, gg=data_paths[0], args=args,
                                                  execute=worker.submit)
            except Exception as exc:  # e.g. an LLM error — keep the session alive
                print(f"canopy: {exc}", file=sys.stderr)
                continue
            if err:
                print(err, file=sys.stderr)
            else:
                sys.stdout.write(out)
                sys.stdout.flush()
            log.say(f"Answered in {gen_s + enh_s + exec_s:.2f}s — llm {gen_s:.2f}s"
                    + (f" · rE2G {enh_s:.2f}s" if enh_s else "")
                    + f" · grove {exec_s:.3f}s")
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
        # The grove is cohort-independent (GENCODE + cCREs). Enhancers are resolved
        # per question from the model's declared COHORT/LAYERS — see _answer.
        _prepare()
        resources_block, preamble, data_paths = _grove_context()
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

    try:  # one-shot: a fresh sandbox per invocation
        out, err, _gen_s, _enh_s, _exec_s = _answer(
            args.question, system_prompt=system_prompt, base=preamble, gg=data_paths[0],
            args=args,
            execute=lambda s: sandbox.run(s, data_paths=data_paths, extra_syspath=[site_dir]),
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
