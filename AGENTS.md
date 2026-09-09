# AGENTS.md

This file provides guidance to any coding agent (Claude Code, Codex, etc.) working
in this repository.

## Project overview

`canopy` is the natural-language interface to [genogrove](https://github.com/genogrove/genogrove).
It is a pure-Python application, **not** the C++ core: given a question, it asks an
LLM to write Python against the [pygenogrove](https://github.com/genogrove/pygenogrove)
bindings, runs that generated code in a sandboxed subprocess, and prints the answer.
There is no compile step for this repo itself.

```
  question ──▶ llm.py ──▶ generated Python ──▶ sandbox.py ──▶ answer
                 ▲                                  │
                 │  pinned datasets + build         │  no network, allowlisted
                 └────────── resources.py ──────────┘  imports, resource caps
```

Dependency direction is one-way — `canopy → pygenogrove` — so `pip install pygenogrove`
never drags in an LLM SDK.

## Key modules

- **`cli.py`** — thin argument-parsing wrapper around the pipeline below.
- **`llm.py`** — Anthropic codegen: builds the prompt from `prompts/system.md` +
  the layer catalogue, sends the question, parses cohort/target declarations and
  the generated Python out of the reply.
- **`sandbox.py`** — **security-critical.** Runs untrusted, LLM-generated code in a
  separate process: network blocked, imports allowlisted, filesystem reads
  restricted to pinned dataset paths, wall-clock/memory/output caps enforced.
  Also owns the warm-grove reuse pattern (`Worker`/`_CANOPY_STATE`) for
  `-i/--interactive` and `serve`.
- **`resources.py`** — the Level 2 reproducibility layer: every dataset is a
  pinned, sha256-verified download from an immutable reference (a commit, never
  a branch/`resolve/main`). Also carries `_GROVE_SCHEMA` and the
  `PYGENOGROVE` build pin.
- **`gff.py`** — GENCODE GFF3 → the universal `Grove` backbone (genes, transcripts,
  exons, splice-chain edges).
- **`layers/`** — one module per data layer (`ccres.py`, `enhancers.py`, `sv.py`),
  each exporting a `LAYER: Layer` descriptor (`_base.py`). `layers/__init__.py`
  collects them into `REGISTRY`, which is the single source of truth the system
  prompt is generated from — adding a layer means adding a module here and
  listing it in `REGISTRY`, nothing else.
- **`prompts/system.md`** — the codegen system prompt: the model's actual API
  contract. This is code, not prose docs — a change here needs the same care as
  a change to `sandbox.py`'s allowlist, and `tests/test_api_surface.py` /
  `tests/test_resources_pins.py` guard it against drifting from the pinned build.
- **`serve.py`** — local web front end over the same pipeline.

## Conventions specific to this project

- **Pin discipline**: a `Resource` in `resources.py` is only ever fetched by
  URL **and** verified by sha256. Artifact URLs must name an immutable
  reference. Never add a resource that resolves against a moving branch/tag.
  If you change what a pinned grove's *content* looks like (payload shape,
  which types are indexed vs. external, new edges), bump `_GROVE_SCHEMA` in the
  same change — that forces every cached copy to be treated as stale and
  re-fetched.
- **Build-time vs. runtime dependencies**: if a transform only needs to run
  once to produce a pinned artifact (e.g. a liftover, a format conversion),
  write it as a one-off script under `tools/` and pin the *derived* output —
  don't wire the transform's own dependency into the package's runtime
  `dependencies` in `pyproject.toml`. See `tools/liftover_pcawg_sv.py` for the
  pattern.
- **Per-question layers** (rE2G enhancer cohorts, PCAWG SV cohorts) attach
  directly to an already-deserialized `Grove` inside the sandbox, **additively**:
  a warm session's grove holds every cohort asked about so far, edges carry the
  cohort (and, for SVs, the sample), and the generated code filters on the
  `COHORTS` / `SV_COHORTS` the host echoes for the current question. Never bake
  cohort-specific structure into the shared pinned artifact. `sv.attach_tracked`
  still returns what it created and `detach` removes exactly that, for callers
  that need a clean working copy.
- **No AI attribution** in commits, PR bodies, or issues — no `Co-Authored-By`
  trailer, no "Generated with Claude Code" footer.
- Don't add abstractions, config, or error handling for scenarios that can't
  happen here — this codebase favors the shortest correct diff over
  speculative generality.

## Build & test

No compiled extension in this repo — `pygenogrove` (the compiled dependency) is
pinned in `pyproject.toml` / `[tool.uv.sources]`.

```bash
uv sync --extra dev                     # full env, including a real pygenogrove build
uv run pytest -q                        # full suite against the real bindings
PYTHONPATH=src uv run --no-project --with pytest pytest -q   # fast skeleton suite, no pygenogrove build
```

CI runs the skeleton suite on every PR; a path-filtered job (`api-surface.yml`)
compiles the pinned `pygenogrove` and runs the bindings-dependent tests when the
pin or the codegen contract changes (see that workflow's own comments for why).

A stray `pytest -q` failure in `test_system_prompt_target_matches_registry`
about a `pygenogrove` version mismatch is a pre-existing, tracked issue, not
something a normal change should fix incidentally — bump the pin deliberately
if you're touching it.

## Where the roadmap lives

`README.md` — install/use instructions, the "what you can ask about" layer
summary, and the reproducibility model — is the user-facing source of truth;
keep it in sync with `layers/__init__.py:REGISTRY` when a layer is added or
its scope changes.
