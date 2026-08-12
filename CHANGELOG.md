# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres
to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Initial project skeleton: `uv` + hatchling packaging, GPL-3.0-or-later license,
  module layout (`cli`, `llm`, `sandbox`, `registry`, `prompts/system.md`), CLI
  smoke tests, and a CI workflow stub.
- **`pygenogrove` API surface in the codegen system prompt** (`prompts/system.md`):
  the bound `Grove` / `BedGrove` / `GffGrove` surface (insert / intersect / flanking),
  the directed-edge graph overlay, file readers, and serialization — with the
  coordinate-convention rules (closed `Interval` vs half-open BED vs 1-based GFF) and
  a worked 2-hop connected-interval example. Documented against `pygenogrove` 0.2.0.
- **Build pinning for Level 2 reproducibility** (`registry.py`): `BuildPin`,
  `verify_pygenogrove_build()` (raises on version drift, returns the C++ engine
  version), and `build_manifest()` for run provenance. `pygenogrove` is pinned to the
  immutable commit `1a9c975` (tag `v0.2.0`) in `pyproject.toml` and mirrored here.
- **Pin-drift guard test** (`tests/test_registry_pins.py`): asserts the `pyproject.toml`
  `==<version>` pin and the `[tool.uv.sources]` `rev` both match `registry.PYGENOGROVE`,
  and that the pin is a full immutable commit SHA — so the Level 2 "all three agree"
  check fails CI on drift instead of relying on manual review. Parses `pyproject.toml`
  with regexes (no `tomllib`) so it runs on the py3.9 floor
  ([#2](https://github.com/genogrove/ask/pull/2)).
- **Out-of-process sandbox** (`sandbox.py`): runs untrusted model-generated Python in an
  isolated subprocess with parent/OS-enforced hard guarantees (stripped env, `setrlimit`
  CPU/memory/no-write/fd caps, whole-session wall-clock kill, byte-capped output) plus
  in-child defense-in-depth (import allowlist with the network/exec primitives scrubbed,
  read-only `open` restricted to registry data roots). The hard boundary is the parent/OS
  layer; an OS-level backend (seccomp/namespaces) is the documented next step for
  adversarial robustness. Covered by 22 isolation tests
  ([#3](https://github.com/genogrove/ask/pull/3)).
- **Baked cCRE grove, per-question enhancers, and `canopy serve`**: the shipped grove is now
  cohort-independent — GENCODE structure plus the ENCODE cCRE registry baked in as
  `regulatory_region` (SO:0005836) nodes, so one `intersect` returns genes and cCREs from a
  single handle. The rE2G enhancer→gene layer moved out of the grove: the model declares
  `COHORT`/`TARGETS` alongside its program, the host grounds the cohort against the ENCODE
  catalog (never silently substituting an unmatched tissue), fetches only the needed links via
  tabix, and injects them as `ENHANCERS`. Adds `canopy serve`, a local web front-end over the
  same generate → execute → render pipeline on a stdlib `http.server`, streaming ndjson
  progress events so the download/build stages are visible. Enhancer records carry
  `ccre_overlap` as a list with per-cCRE shared base counts, because most rE2G windows span
  several cCREs and about a third span differing classes
  ([#5](https://github.com/genogrove/canopy/pull/5)).

### Changed
- **Upgraded `pygenogrove` to v0.4.0** (the universal JSON Grove redesign): re-pinned to the
  immutable commit `d6c75b9` (tag `v0.4.0`) across `pyproject.toml` + `registry.PYGENOGROVE`,
  and rewrote the `prompts/system.md` codegen contract for the new surface — `GenomicCoordinate`
  is the one key type (`Interval` removed; strand-aware, `'*'` wildcard query), `Grove` is the
  universal `grove<genomic_coordinate, json>` storing arbitrary JSON payloads, plus strand-aware
  `intersect`/`flanking` (predicate overload), `remove_key`/`compact`/counts, and typed
  `BedGrove`/`GffGrove` for C++ interop ([#4](https://github.com/genogrove/ask/pull/4)).
- **Harmonized the import package with the distribution name**: `canopy` → `genogrove_canopy`,
  matching the already-declared `genogrove-canopy` dist. PyPI `canopy` is taken by an actively
  released project (9.10) that installs its own top-level `canopy/` directory, so co-installing
  the two silently clobbered one or the other in `site-packages`. Flat rather than a
  `genogrove.canopy` namespace package, since the bindings are a compiled top-level extension
  module that cannot join a namespace without being rebuilt. The product name, the `canopy`
  console script, and the paper's `genogrove canopy <question>` surface are unchanged; existing
  environments need `uv sync` to regenerate the console script
  ([#6](https://github.com/genogrove/canopy/pull/6)).

### Fixed
- **One unified grove, pinned from Hugging Face**: the shipped artifact is now the GENCODE
  backbone with the 2,348,854 ENCODE cCREs already built into it, pinned at an immutable
  `genogrove/canopy` commit. This fixes a cold install that was broken for everyone but the
  author — the cCRE resource carried `pending-upload.invalid` URLs and worked only from a
  pre-seeded cache, so a fresh `--init` died on DNS while baking — and replaces a pinned grove
  built under the old payload model, against which every generated query raised `KeyError` on
  the first chain filter. The local bake is gone (`ensure_baked_grove`, `_baked_grove_gg`,
  `_BAKED_SCHEMA`), along with the deserialize-then-insert workaround for pygenogrove#68.
  `_GROVE_SCHEMA` 2 → 3 so cached copies of the old artifact are re-fetched rather than
  silently reused. The pin guard now also covers the data pins: full-length checksums, no
  movable refs, and no placeholder URLs on anything resolved at runtime
  ([#9](https://github.com/genogrove/canopy/issues/9), [#11](https://github.com/genogrove/canopy/pull/11)).
- **`prompts/system.md` declared the wrong pinned build**: the header said "Current target:
  pygenogrove 0.6.2" while the pin was 0.7.4. Because `system.md` *is* the system prompt, that
  line shipped to the model on every request, telling it the documented API surface targeted a
  build three releases behind the one its generated code runs against. `test_resources_pins.py`
  now asserts the header version against `resources.PYGENOGROVE`, closing the one leg of the
  Level 2 "all agree" invariant that had no guard
  ([#7](https://github.com/genogrove/canopy/issues/7), [#8](https://github.com/genogrove/canopy/pull/8)).
