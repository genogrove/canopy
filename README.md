# canopy

**Ask plain-English questions over connected genomic intervals.**

`canopy` is the natural-language interface to [genogrove](https://github.com/genogrove/genogrove).
You ask a question; it writes Python against the [pygenogrove](https://github.com/genogrove/pygenogrove)
bindings, runs that code in a sandbox, and prints the answer.

```console
$ canopy "What enhancers regulate AR in prostate cancer?"
```

```
AR enhancers in prostate cancer (LNCaP): 22 link(s), 16 distal:
strand=.  type=enhancer  n=1  cohort=EFO:0005726  target=AR
chrom  start     end       name                class       score    ccre_overlap
chrX   67544302  67545353  enh:promoter->AR    promoter    0.99999  pELS:332 PLS:277 PLS:162
chrX   67545965  67547497  enh:genic->AR       genic       0.99888  pELS:329 pELS:318 dELS:214
chrX   67551264  67551763  enh:genic->AR       genic       0.6421   -
chrX   67538334  67538833  enh:intergenic->AR  intergenic  0.47963  TF:168
…
```

Values identical on every row are stated once above the table. `-` means *no* overlapping
cCRE — a real finding for about 2% of links, not missing data. `--format json` keeps the full
structure, including the cCRE accessions.

## Why

genogrove stores annotations as a *connected* structure: intervals in per-chromosome
B+ trees with a directed graph overlay linking related keys — exon→transcript,
breakpoint→mate, enhancer→gene. Questions that would otherwise need a brittle
`intersect | awk | sort | join` pipeline become one traversal.

`canopy` makes that structure reachable without writing code. There is no fixed query
language: the model targets the bindings directly, so anything the bindings can express
is askable.

## Install

`pygenogrove` is a C++/htslib extension built from source, so you need a compiler, CMake
and htslib. The project uses [`uv`](https://docs.astral.sh/uv/).

```console
$ brew install uv htslib cmake          # macOS; on Linux use your package manager
$ git clone https://github.com/genogrove/canopy && cd canopy
$ env CMAKE_PREFIX_PATH=/opt/homebrew \
      CMAKE_ARGS="-DCMAKE_PREFIX_PATH=/opt/homebrew/opt/htslib" uv sync
```

Then fetch the data once (a pinned ~109 MB grove) so the first question is instant:

```console
$ uv run canopy --init
```

## Use

Set `ANTHROPIC_API_KEY`, then ask. (`uv run` is only needed until you activate the
environment — `source .venv/bin/activate`, after which `canopy` works on its own.)

```console
$ uv run canopy "Which gene contains the variant at chr7:55,191,822?"
$ uv run canopy --show-code "Which transcripts share an exon with chr7:55,191,822?"
$ uv run canopy --format json "List the exons of EGFR"
```

| Flag | |
|---|---|
| `--format text\|bed\|tsv\|json` | Output format. `text` is an aligned table; the rest are machine-readable. |
| `--show-code` | Print the generated Python before running it. |
| `-i, --interactive` | Keep the grove open across questions — the ~200 ms open is paid once, then queries are sub-millisecond. |
| `--cohort NAME` | Pick the ENCODE-rE2G biosample for enhancer questions (`--list-cohorts` to browse 369 of them). |
| `--init` | Download the grove now and exit. |

There is also a local web front end over the same pipeline:

```console
$ uv run canopy serve
```

## What you can ask about

One grove, queried through one handle, holds three layers:

- **Gene structure** — GENCODE v50: genes, transcripts and exons, with `contains` and
  `first_exon`/`next` splice-chain edges, and CDS ranges on each exon.
- **Candidate regulatory elements** — the ENCODE cCRE registry (V4, 2,348,854 elements),
  typed as Sequence Ontology `regulatory_region` with the evidence-based class (`PLS`,
  `pELS`, `dELS`, …) in the payload. A single `intersect` returns genes *and* cCREs.
- **Enhancer→gene links** — ENCODE-rE2G predictions across 369 biosamples. These are
  cohort-specific, so they are resolved per question rather than shipped in the grove.

Enhancer answers carry a `ccre_overlap` list rather than a single class: most rE2G
windows span several cCREs, and about a third span cCREs of differing classes, so
collapsing them to one would assert a link the data does not make.

## Reproducibility

canopy targets **Level 2** reproducibility — every dataset and every library build is
pinned:

- Datasets are pinned by URL **and** sha256, verified on download. Only names in the
  curated catalog are ever fetched; open-web resource discovery is out of scope.
- Artifact URLs must name an immutable reference (a commit, not a branch), so the bytes
  behind a pin cannot change underneath the checksum.
- The `pygenogrove` build is pinned to a commit SHA, mirrored in `pyproject.toml` and in
  the codegen prompt, with a test asserting all of them agree.

Given the same question and the same catalog, a run reproduces.

## Configuration

| Variable | |
|---|---|
| `ANTHROPIC_API_KEY` | Required for asking questions. Not needed for the data layer. |
| `GENOGROVE_CANOPY_CACHE` | Where datasets are cached. Defaults to `~/.cache/genogrove-canopy`. Point it at a temporary directory to exercise a cold fetch without touching the real cache. |

## How it works

```
  question ──▶ llm.py ──▶ generated Python ──▶ sandbox.py ──▶ answer
                 ▲                                  │
                 │  pinned datasets + build         │  no network, allowlisted
                 └────────── resources.py ──────────┘  imports, resource caps
```

The generated code is untrusted, so it never runs in-process: `sandbox.py` executes it
in a separate process with the network blocked, imports allowlisted, filesystem reads
restricted to pinned dataset paths, and wall-clock, memory and output caps enforced.

Dependency direction is one-way — `canopy → pygenogrove` — so `pip install pygenogrove`
never drags in an LLM SDK.

## The `genogrove canopy` surface

In the genogrove paper and docs the command is written `genogrove canopy <question>`.
That is a thin alias over the `canopy` console script this package installs; the
application layer lives here, separate from the core C++ CLI, because it has a different
release cadence and audience.

## Development

```console
$ uv run --extra dev pytest -q          # full suite against the real bindings
```

CI runs a fast skeleton suite on every PR, plus a job that compiles the pinned
`pygenogrove` and runs the bindings-dependent tests when the pin or the codegen contract
changes.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
