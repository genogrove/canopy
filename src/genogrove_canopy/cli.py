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
from pathlib import Path

from genogrove_canopy import __version__, llm, log, resources, sandbox

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
QUERY_SURFACE = ("intersect", "flanking", "get_neighbors", "get_edges", "get_neighbors_if")


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


def _resolve_cohorts(specs):
    """Map ``--cohort`` specs to ``{cohort name: [accessions]}`` via ``re2g_cohorts``.

    Each spec matches a cohort by exact ontology id or case-insensitive substring of its
    name; the most-replicated match wins (the catalog is sorted by replicate count). Raises
    ``SystemExit`` with a pointer to ``--list-cohorts`` on no match.
    """
    catalog = resources.re2g_cohorts()
    chosen = {}
    for spec in specs:
        s = spec.strip().lower()
        match = next((c for c in catalog if c["ontology_id"].lower() == s), None) or \
            next((c for c in catalog if s in c["name"].lower()), None)
        if match is None:
            raise SystemExit(f"canopy: no cohort matches {spec!r} — see --list-cohorts")
        chosen[match["name"]] = match["accessions"]
    return chosen


def _list_cohorts() -> None:
    """Print the available rE2G cohorts (most-replicated first) to stdout."""
    print(f"{'ontology id':16}  {'reps':>4}  {'type':16}  name")
    for c in resources.re2g_cohorts():
        print(f"{c['ontology_id']:16}  {c['n_replicates']:>4}  {c['type']:16}  {c['name']}")


def _grove_context():
    """Resolve the shipped grove to ``(resources_block, code_preamble, data_paths)``.

    The grove is the GENCODE backbone with the **Tier-1 static layers already in it** — currently
    the ENCODE cCRE registry, built into the pinned artifact rather than baked on first run (see
    ``resources.ensure_all_grove``). So a `intersect`
    returns genes *and* cCREs from one handle, ``GENCODE_HUMAN``, opened lazily. The enhancer layer
    is **not** in the grove (it is dynamic/cohort-specific): the host resolves the model's declared
    ``COHORT``/``TARGETS`` through the tabix index and injects only the needed enhancers as the
    ``ENHANCERS`` variable (defaulted to ``[]`` in the preamble so the code never ``NameError``s).
    """
    from genogrove_canopy import layers

    var = "GENCODE_HUMAN"
    gg = str(resources.ensure_all_grove(_BASE))
    block = resources_block(
        var, resources.RESOURCES[_BASE].description, layers.catalogue_block(["ccre"])
    )
    preamble = f"{var} = {json.dumps(gg)}\nENHANCERS = []\n"
    return block, preamble, [gg]


def resources_block(var: str, description: str, layers_block: str) -> str:
    """Render the "Available resources" text the model reads, for grove handle ``var``.

    Pure so it can be tested without a grove: ``_grove_context`` resolves the ~109 MB artifact,
    which a unit test has no business downloading. That matters here because this is the only place
    ``QUERY_SURFACE`` reaches the model — a test asserting those names appear in the *rendered*
    block is checking the real contract, not that a token exists somewhere in a source file.
    """
    return (
        f"- `{var}` (str): path to the shipped grove "
        f"({description}) — gene/transcript/exon structure **plus the "
        f"ENCODE cCRE nodes, in the same grove**. Open it lazily with `g = pg.GroveView.open({var})`. A "
        f"**located** query (a variant at chr7:55191822) reads just that locus; a **genome-wide / "
        f"gene-name** query works from the same handle. Query-only: "
        f"{', '.join(f'`{m}`' for m in QUERY_SURFACE)}.\n"
        f"  Node layers in the grove — returned by `intersect` alongside genes, filter on `type`:\n"
        f"  {layers_block}\n"
        f"- `ENHANCERS` (list): the ENCODE-rE2G enhancer→gene links for the `COHORT`/`TARGETS` you "
        f"declare (see \"Enhancers\"). Empty `[]` unless the question is about enhancers/regulation."
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


def _resolve_query_cohorts(args, cohort_hint):
    """Grounded cohort resolution for a query that needs enhancers, precedence-ordered:
    ``--cohort`` override → the model's declared ``COHORT`` (matched against the catalog; no
    match → none, so a wrong tissue is never silently substituted) → the default cohort.
    Returns ``({name: accessions}, note)`` where ``note`` is a stderr line or ``None``."""
    if args.cohort:
        return _resolve_cohorts(args.cohort), None
    if cohort_hint:
        try:
            return _resolve_cohorts([cohort_hint]), None
        except SystemExit:  # the model named a tissue with no catalog match — don't substitute
            return {}, f"no ENCODE cohort matched {cohort_hint!r} — no enhancers loaded (see --list-cohorts)"
    return _resolve_cohorts([DEFAULT_COHORT]), "default"


def _prepare() -> None:
    """First-run notice if the shipped grove isn't cached yet: a one-time ~109 MB download of the
    pinned unified `.gg` (gene structure + cCREs). No local build."""
    if not resources._all_grove_gg(_BASE).exists():
        log.say(f"Fetching the grove ({resources.RESOURCES[_BASE].grove_contents}) — "
                "first run only, a pinned ~109 MB .gg")


def _cohort_ids(cohorts) -> list:
    """Ontology ids for the selected cohort **names** — the enhancer index is keyed by id."""
    if not cohorts:
        return []
    name2id = {c["name"]: c["ontology_id"] for c in resources.re2g_cohorts()}
    return [name2id[n] for n in cohorts if n in name2id]


def _answer(question, *, system_prompt, preamble, args, execute):
    """Translate one question to code, run it via ``execute(script)``, and render.

    ``execute`` is a ``script -> SandboxResult`` callable (``sandbox.run`` for one-shot,
    ``Worker.submit`` for interactive). Returns ``(rendered_stdout, error_msg, gen_s, exec_s)``
    — exactly one of stdout/error is non-empty; the two times split code-gen from execution.

    The enhancer layer is resolved **per question**: the model declares ``COHORT``/``TARGETS``,
    the host grounds the cohort (``--cohort`` overrides), fetches only those enhancers via the
    tabix index, and injects them as ``ENHANCERS`` — no whole-cohort grove augment.
    """
    t0 = time.perf_counter()
    cohort_hint, targets, code = llm.generate_query(question, system_prompt, model=args.model)
    gen_s = time.perf_counter() - t0
    if args.show_code:
        print("# --- generated code ---", file=sys.stderr)
        print(code, file=sys.stderr)
    enh_pre = ""
    if targets:  # an enhancer/regulation question — resolve the cohort and fetch its enhancers
        from genogrove_canopy.layers import enhancers
        cohorts, note = _resolve_query_cohorts(args, cohort_hint)
        cohort_ids = _cohort_ids(cohorts)
        records = enhancers.fetch_for_targets(targets, cohort_ids) if cohort_ids else []
        if records:
            enh_pre = enhancers.preamble(records)
            src = " (default — name a tissue or pass --cohort)" if note == "default" else ""
            log.say(f"Enhancers: {len(records)} links from cohort(s) "
                    f"{'; '.join(cohorts)}{src}")
        elif note and note != "default":
            log.say(note)
    # JSONL is the output contract, so guarantee `json` is importable even if the
    # generated code forgets the import (it's already in the allowlist).
    t1 = time.perf_counter()
    result = execute("import json\n" + preamble + enh_pre + code)
    exec_s = time.perf_counter() - t1
    if result.returncode != 0 or result.timed_out:
        return "", (result.stderr.strip() or "(the generated code failed with no output)"), gen_s, exec_s
    rendered = _render(result.stdout, args.format)
    if not rendered.strip():
        return "", "(the generated code produced no output)", gen_s, exec_s
    return rendered, "", gen_s, exec_s


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
                out, err, gen_s, exec_s = _answer(question, system_prompt=system_prompt,
                                                  preamble=preamble, args=args, execute=worker.submit)
            except Exception as exc:  # e.g. an LLM error — keep the session alive
                print(f"canopy: {exc}", file=sys.stderr)
                continue
            if err:
                print(err, file=sys.stderr)
            else:
                sys.stdout.write(out)
                sys.stdout.flush()
            log.say(f"Answered in {gen_s + exec_s:.2f}s — llm {gen_s:.2f}s · grove {exec_s:.3f}s")
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
        # per question from the model's declared COHORT/TARGETS — see _answer.
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
        out, err, _gen_s, _exec_s = _answer(
            args.question, system_prompt=system_prompt, preamble=preamble, args=args,
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
