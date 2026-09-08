# SPDX-License-Identifier: GPL-3.0-or-later
"""Curated resource catalog — the Level 2 reproducibility layer.

A run is reproducible when the question, the resolved datasets, and the library
builds are all pinned. This module is the single source of truth for those pins:

* **Datasets** — named genomic resources with a pinned URL and checksum, resolved
  to a local path on demand (with checksum verification).
* **Builds** — the exact ``pygenogrove`` / ``genogrove`` versions a run was made
  against, recorded so results can be regenerated.

Open-web resource discovery (Level 3) is intentionally out of scope. ``resolve``
takes a curated *name*, never a URL — so the only data ever fetched is what
``RESOURCES`` explicitly defines.
"""

from __future__ import annotations

import hashlib
import os
import shlex
import shutil
import subprocess
import tempfile
import urllib.request
from collections.abc import Iterable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from urllib.parse import urlsplit


# --------------------------------------------------------------------------- #
# Builds — the pinned library versions a run is made against.
#
# pygenogrove is resolved from a pinned git commit (see [tool.uv.sources] in
# pyproject.toml). The pin below mirrors that commit so a run can record, and
# verify against, the exact build it was made with.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class BuildPin:
    """An exact, immutable library build a run is reproducible against."""

    name: str
    version: str  # expected package version (``pygenogrove.__version__``)
    git_rev: str  # immutable commit the build is pinned to
    git_tag: str = ""  # human-readable tag at that commit, if any


# Keep in lockstep with [tool.uv.sources] `rev` in pyproject.toml, and with the "Current target:"
# line in prompts/system.md — that header ships inside the system prompt, so a stale version tells
# the model its documented API surface targets a build its generated code won't run against. All
# four are asserted against each other by tests/test_resources_pins.py; bump them in one commit.
PYGENOGROVE = BuildPin(
    name="pygenogrove",
    version="0.9.0",
    git_rev="3846c5de0084b08e91a724eee5021656d6148290",
    git_tag="v0.9.0",
)


def verify_pygenogrove_build() -> str:
    """Check that the installed ``pygenogrove`` matches the pinned build.

    Imports lazily so this module is usable without ``pygenogrove`` installed
    (e.g. in the skeleton test env). Returns the underlying C++ engine version
    (``pygenogrove.__genogrove_version__``) so a run can record it. Raises
    ``RuntimeError`` on version drift from the pin.
    """
    import pygenogrove

    installed = getattr(pygenogrove, "__version__", None)
    if installed != PYGENOGROVE.version:
        raise RuntimeError(
            f"pygenogrove build drift: pinned {PYGENOGROVE.version} "
            f"(rev {PYGENOGROVE.git_rev}), installed {installed!r}. "
            "Run `uv sync` to match the pinned build."
        )
    return str(getattr(pygenogrove, "__genogrove_version__", ""))


def build_manifest() -> dict[str, str]:
    """Provenance record of the build a run was made against.

    Combines the static pin with the engine version observed at runtime, for
    embedding in a run's output so results can be regenerated.
    """
    return {
        "pygenogrove_version": PYGENOGROVE.version,
        "pygenogrove_git_rev": PYGENOGROVE.git_rev,
        "pygenogrove_git_tag": PYGENOGROVE.git_tag,
        "genogrove_engine_version": verify_pygenogrove_build(),
    }


# --------------------------------------------------------------------------- #
# Datasets — pinned genomic resources (URL + checksum), resolved to a local path.
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class Resource:
    """A pinned genomic dataset in the curated catalog.

    ``url`` + ``sha256`` pin the annotation file. For region access it must be
    **bgzip-compressed and tabix-indexed**; when the ``.tbi`` is hosted too,
    ``index_url`` + ``index_sha256`` pin it and ``resolve``/``indexed_path`` fetch
    the pair (no local indexing). ``filename`` overrides the local name when it
    can't be derived from the URL (e.g. Zenodo's ``…/files/<name>/content``).
    """

    name: str
    url: str
    sha256: str
    description: str = ""
    filename: str = ""
    index_url: str = ""
    index_sha256: str = ""
    # Prebuilt whole-genome grove (.gg). When set, the genome-wide path downloads this
    # instead of building locally (minutes) — GroveView reads it lazily. Rebuild + re-pin
    # if the .gg format or the source annotation changes.
    grove_url: str = ""
    grove_sha256: str = ""
    # Layers the pinned grove carries *beyond* the annotation itself, named for humans. A local
    # build reads only `url`, so it cannot reproduce these — declaring them makes
    # `ensure_all_grove` refuse to substitute a structurally different grove instead of silently
    # returning one whose extra layers are simply absent (queries would come back empty, not error).
    grove_layers: tuple[str, ...] = ()
    # Short human label for progress messages — what is *in* the grove, not the catalog key.
    # `gencode.human grove` reads as GENCODE-only and understates an artifact that also carries
    # 2.35M cCREs; `grove_layers` above is the machine-checked list and is too long for one line.
    grove_contents: str = ""


# Curated dataset catalog. Each entry pins an *immutable* release (an explicit
# version, never a "latest"/"current" symlink) by URL + sha256. Only names listed
# here are ever fetched.
RESOURCES: dict[str, Resource] = {
    "gencode.human": Resource(
        name="gencode.human",
        # Coordinate-sorted, bgzip+tabix build of GENCODE v50 (Zenodo 21123308),
        # derived from GENCODE v50 (upstream sha 2aaf245c…875899a) — see the record.
        url="https://zenodo.org/api/records/21123308/files/gencode.v50.annotation.sorted.gff3.gz/content",
        sha256="2a87d3a39f9e3be6f0c49359724223ba5e0a094f2fc059b2655635888bb223f5",
        filename="gencode.v50.annotation.sorted.gff3.gz",
        index_url="https://zenodo.org/api/records/21123308/files/gencode.v50.annotation.sorted.gff3.gz.tbi/content",
        index_sha256="52020642c93f01c24488d98b446d705a655d31ea39339fad36cced3b9cc9480a",
        # Prebuilt **unified** grove (105 MB): the GENCODE backbone — indexed types are exactly
        # `{gene, regulatory_region}`, verified directly (transcript is EXTERNAL, like exon:
        # reach one via its gene's `contains` edge, never via intersect(); coding-structure
        # annotations, including the selenocysteine-readthrough analogue of `stop_codon`, are
        # dropped like every other UTR/codon type, never indexed) — PLUS the 2,348,854 ENCODE
        # cCREs as `type:"regulatory_region"`, `source:"ENCODE-SCREEN"` nodes. No SV-layer-
        # specific structure baked in — `layers/sv.py` computes its own segment/gap coverage on
        # demand, per session, into an already-deserialized working copy instead (verified safe:
        # inserting into a deserialized grove does not corrupt it). 2,427,587 indexed keys /
        # 4,204,002 keys total (indexed + external), built in one pass via `gff.load_gff` +
        # `layers.ccres.attach`.
        #
        # Rebuilt (content unchanged, only the on-disk format changed) under the actually-pinned
        # pygenogrove 0.9.0 — the previous pin was built under a stale local 0.7.4 install despite
        # pyproject.toml already naming 0.9.0, so it silently shipped an old serialization format
        # (0.2) that 0.9.0's `GroveView.open`/`Grove.deserialize` reject outright ("bad magic (not
        # a format 0.3 grove stream)"). Round-tripped and EGFR-region smoke-checked against the new
        # build before re-pinning.
        #
        # Because the layer ships inside the grove, nothing resolves `encode.ccre.v4` at query time
        # and there is no local bake — the old deserialize→insert→reserialize step (and the
        # `_BAKED_SCHEMA` that versioned it) is gone.
        #
        # The URL pins an immutable HF **commit**, not `resolve/main`: a branch ref is movable and
        # would defeat the Level 2 guarantee. Verified end-to-end through `_download` (302 → 200,
        # 104,994,275 bytes, sha256 match).
        grove_url=(
            "https://huggingface.co/datasets/genogrove/canopy/resolve/"
            "826d9f99599b25a81f70b30f592d8ba89d3d0c20"
            "/groves/gencode.v50+ccre.v4.grove-model.gg"
        ),
        grove_sha256="1b3d12b23c8d218ed2c83c38f3ba1e5c7dbad2cbd1d668874e6c4725b0ed51a9",
        grove_layers=("ENCODE cCRE registry (V4, GRCh38)",),
        grove_contents="GENCODE human v50 + ENCODE cCREs V4",
        description="GENCODE v50 comprehensive gene annotation, GRCh38 (GFF3, sorted + bgzip + tabix).",
    ),
    "encode.ccre.v4": Resource(
        name="encode.ccre.v4",
        # ENCODE Registry of cCREs V4 (GRCh38). Nature 2026, doi:10.1038/s41586-025-09909-9.
        #
        # **Build-time input — not resolved on a user's machine.** The cCREs ship inside the
        # unified grove above, so the query path reads them from there; nothing in `src/` resolves
        # this entry (`layers.ccres.all_records`/`in_region` are host-side helpers with no
        # production caller). It is pinned anyway so *rebuilding* the unified grove is reproducible
        # rather than depending on one machine's cache — which is what it did until now, when both
        # URLs were `pending-upload.invalid` placeholders and these 32 MB existed nowhere else.
        #
        # Uploaded from the content-addressed cache, so the bytes are the ones the sha256s below
        # already pinned; both were re-fetched through `_download` and verified.
        url=(
            "https://huggingface.co/datasets/genogrove/canopy/resolve/"
            "23fa3a6420cba225c7ba8821d11530a404937655"
            "/ccres/GRCh38-cCREs.v4.bed.gz"
        ),
        sha256="9f33d157de568afffedc0b3bebd0b5aaa350e341cb18d5009d2fee6f4a8cee0d",
        filename="GRCh38-cCREs.v4.bed.gz",
        index_url=(
            "https://huggingface.co/datasets/genogrove/canopy/resolve/"
            "23fa3a6420cba225c7ba8821d11530a404937655"
            "/ccres/GRCh38-cCREs.v4.bed.gz.tbi"
        ),
        index_sha256="485e642dd8f0fb97ff693157a54ef47e881cb31ee8dcfca54bea73bfb64721cd",
        description="ENCODE Registry of cCREs V4, GRCh38 (2,348,854 elements; Nature 2026, "
                    "doi:10.1038/s41586-025-09909-9).",
    ),
    "pcawg.sv.icgc": Resource(
        name="pcawg.sv.icgc",
        # PCAWG consensus structural-variant calls (v1.6), ICGC portion. Open access, no DACO —
        # confirmed against the source's own README (open/controlled split lives at the file
        # level, this one is public). One BEDPE per sample inside the tarball
        # (`icgc/open/<aliquot_id>.pcawg_consensus_1.6.161116.somatic.sv.bedpe.gz`), 1,926 samples.
        #
        # PCAWG's own release is hg19 only — no official hg38 SV release exists — so this is a
        # *derived* artifact: lifted to GRCh38 by `tools/liftover_pcawg_sv.py` from the pinned
        # `pcawg.sv.icgc.hg19` original (build-time step, not resolved on a user's machine; same
        # shape as `encode.ccre.v4`) and re-hosted, since there's no upstream hg38 original. 179,905/179,973 SVs lifted (0.04%
        # dropped — didn't lift cleanly). Chrom columns normalized to "chr"-prefixed to match the
        # GENCODE backbone. Re-lifted once with the strand of reverse-mapped breakends corrected
        # (778 records); the dataset's `layers/genomic/` layout starts with this commit.
        url=(
            "https://huggingface.co/datasets/genogrove/canopy/resolve/"
            "b3abfea82ad29de7dbadda47baeec5a99c4bce77"
            "/layers/genomic/pcawg-sv/final_consensus_sv_bedpe_passonly.icgc.public.hg38.tgz"
        ),
        sha256="b06722c0ed2f7ad22085d2a3630b3711516b7d071313dd4a25b66fd63e1cacd6",
        filename="final_consensus_sv_bedpe_passonly.icgc.public.hg38.tgz",
        description="PCAWG consensus SV calls (v1.6), ICGC portion, 1,926 samples, open access, "
                    "lifted to GRCh38.",
    ),
    "pcawg.sv.tcga": Resource(
        name="pcawg.sv.tcga",
        # Same release, TCGA portion — 822 samples, same open/lifted-and-re-hosted shape.
        # 129,215/129,273 SVs lifted (0.04% dropped).
        url=(
            "https://huggingface.co/datasets/genogrove/canopy/resolve/"
            "b3abfea82ad29de7dbadda47baeec5a99c4bce77"
            "/layers/genomic/pcawg-sv/final_consensus_sv_bedpe_passonly.tcga.public.hg38.tgz"
        ),
        sha256="d7c6fb040b518b8f590a685c7da1de372d4142b923b5da2da06ea8091ed8a5b6",
        filename="final_consensus_sv_bedpe_passonly.tcga.public.hg38.tgz",
        description="PCAWG consensus SV calls (v1.6), TCGA portion, 822 samples, open access, "
                    "lifted to GRCh38.",
    ),
    # The hg19 ORIGINALS the two entries above were lifted from — build inputs for
    # `tools/liftover_pcawg_sv.py`, never read at query time. Hosted on the ICGC 25K open
    # bucket (S3-compatible, public read); the URL carries no version, so immutability rests on
    # the sha256 here: a changed upstream file fails verification instead of lifting silently.
    "pcawg.sv.icgc.hg19": Resource(
        name="pcawg.sv.icgc.hg19",
        url="https://object.genomeinformatics.org/icgc25k-open/PCAWG/consensus_sv/"
            "final_consensus_sv_bedpe_passonly.icgc.public.tgz",
        sha256="8aff040b21a3680629a364e3cfc57f42a235cc8d8aa3e1a14245f5e2c9d01079",
        filename="final_consensus_sv_bedpe_passonly.icgc.public.tgz",
        description="PCAWG consensus SV calls (v1.6), ICGC portion, hg19 original (liftover input).",
    ),
    "pcawg.sv.tcga.hg19": Resource(
        name="pcawg.sv.tcga.hg19",
        url="https://object.genomeinformatics.org/icgc25k-open/PCAWG/consensus_sv/"
            "final_consensus_sv_bedpe_passonly.tcga.public.tgz",
        sha256="3be5502baa2726142137ac83c234623ca215fb356b8645c63933ee516e5c94f6",
        filename="final_consensus_sv_bedpe_passonly.tcga.public.tgz",
        description="PCAWG consensus SV calls (v1.6), TCGA portion, hg19 original (liftover input).",
    ),
}


# Content-addressed cache: <CACHE>/<sha256>/<filename>. A file only lands here
# after its checksum is verified, so a cache hit needs no re-verification.
#
# `GENOGROVE_CANOPY_CACHE` overrides the location, so a cold run can be exercised against a
# throwaway directory instead of deleting the real cache:
#
#     GENOGROVE_CANOPY_CACHE=$(mktemp -d) canopy --init
#
# That is what makes the cold-fetch path testable in CI, which has no cache to clear and must not
# clear a shared one. Read once at import; tests that need to redirect it mid-process monkeypatch
# `resources._CACHE` directly.
_CACHE = Path(
    os.environ.get("GENOGROVE_CANOPY_CACHE") or Path.home() / ".cache" / "genogrove-canopy"
)


def resolve(name: str) -> Path:
    """Resolve a curated resource to a verified local file path.

    ``name`` is a catalog key (``KeyError`` if not curated) — never a URL, so the
    only data ever fetched is what ``RESOURCES`` defines. On a cache miss, streams
    the pinned URL to a temp file while hashing, and only commits it to the cache
    if the sha256 matches. A mismatch is a hard failure: the partial download is
    discarded and nothing is cached.
    """
    res = RESOURCES[name]
    fname = res.filename or Path(urlsplit(res.url).path).name
    return _download(res.url, res.sha256, _CACHE / res.sha256 / fname)


def _download(url: str, sha256: str, dest: Path, label: str = "") -> Path:
    """Stream ``url`` to ``dest`` (cache hit = no-op), verifying its sha256.

    A mismatch is a hard failure: the partial download is discarded, nothing is
    committed. The commit is an atomic rename within ``dest``'s directory. An empty
    ``sha256`` skips verification — for resolve-on-demand files (ENCODE-rE2G) that are
    not pinned until a run freezes them; curated ``RESOURCES`` always pass a checksum.

    ``label`` turns on a progress counter; without one the download is silent, which is what
    library callers and tests want.
    """
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    tmp = tempfile.NamedTemporaryFile(dir=dest.parent, delete=False)
    tmp_path = Path(tmp.name)
    try:
        with tmp, urllib.request.urlopen(url) as resp:  # noqa: S310 — pinned catalog URL
            progress = None
            if label:
                from genogrove_canopy.log import Progress

                total = resp.headers.get("Content-Length")
                progress = Progress(label, int(total) if total and total.isdigit() else None)
            try:
                for chunk in iter(lambda: resp.read(1 << 20), b""):
                    digest.update(chunk)
                    tmp.write(chunk)
                    if progress:
                        progress.advance(len(chunk))
            except BaseException:
                if progress:
                    progress.abort()  # else the error prints onto the half-drawn progress line
                raise
            if progress:
                progress.finish()
        if sha256 and digest.hexdigest() != sha256:
            raise ValueError(
                f"checksum mismatch for {url!r}: expected {sha256}, got {digest.hexdigest()}"
            )
        tmp_path.replace(dest)
    except BaseException:
        tmp_path.unlink(missing_ok=True)
        raise
    return dest


def data_roots(names: Iterable[str]) -> list[str]:
    """Resolve ``names`` to local file paths for the sandbox's read-only roots."""
    return [str(resolve(n)) for n in names]


# --------------------------------------------------------------------------- #
# Region access — bgzip + tabix so a query reads only its locus (see genogrove_canopy.gff.
# build_grove). GENCODE ships plain gzip, so we recompress + index once.
# --------------------------------------------------------------------------- #


def indexed_path(name: str) -> Path:
    """A bgzip-compressed, coordinate-sorted, tabix-indexed GFF for ``name``.

    If the resource pins a hosted index (``index_url``), download the annotation +
    its ``.tbi`` (no local work). Otherwise fall back to building the index locally
    from the plain-gzip download with htslib's ``bgzip``/``tabix`` — a one-time
    ~minutes step; ``RuntimeError`` if those tools are missing.

    Region reads (``pg.GffReader(path, region=...)``) need the ``.tbi`` next to the
    returned ``.gff3.gz``; both paths are placed accordingly.
    """
    res = RESOURCES[name]
    if res.index_url:  # hosted pair — download, don't build
        gff = resolve(name)  # the sorted-bgzip annotation
        _download(res.index_url, res.index_sha256, gff.with_name(gff.name + ".tbi"))
        return gff

    src = resolve(name)  # plain-gzip download
    out = src.with_name("indexed.gff3.gz")
    tbi = out.with_name(out.name + ".tbi")
    if out.exists() and tbi.exists():
        return out
    for tool in ("bgzip", "tabix"):
        if shutil.which(tool) is None:
            raise RuntimeError(
                f"{tool!r} not found — install htslib for region access "
                "(e.g. `brew install htslib` / `apt install tabix`)"
            )
    tmp = out.with_name("indexed.tmp.gff3.gz")
    q = shlex.quote(str(src))
    # Header ('#') lines first, then data sorted by (seqid, start) as tabix requires. A prepended
    # rank column breaks same-start ties into gene -> transcript -> children; without it sort falls
    # back to comparing whole lines, so a gene and its exon sharing a start come out alphabetically
    # (exon, gene, transcript). Harmless for tabix and for the two-pass loader in genogrove_canopy.gff, but
    # the file should read in nesting order. `cut` drops the rank again before bgzip.
    rank = r"""awk -F'\t' -v OFS='\t' '{r=($3=="gene")?1:($3=="transcript")?2:3; print r, $0}'"""
    pipeline = (
        f"{{ gzip -dc {q} | grep '^#' ; "
        f"gzip -dc {q} | grep -v '^#' | {rank} | sort -k2,2 -k5,5n -k1,1n | cut -f2- ; }} "
        f"| bgzip -c > {shlex.quote(str(tmp))}"
    )
    subprocess.run(pipeline, shell=True, check=True)  # noqa: S602 — our own quoted paths
    subprocess.run(["tabix", "-p", "gff", str(tmp)], check=True)
    Path(str(tmp) + ".tbi").replace(tbi)
    tmp.replace(out)  # commit the index (both parts now in place)
    return out


def is_indexed(name: str) -> bool:
    """True if ``name``'s indexed GFF + `.tbi` are already local (no download/build)."""
    res = RESOURCES[name]
    if res.index_url:  # hosted pair
        fname = res.filename or Path(urlsplit(res.url).path).name
        gff = _CACHE / res.sha256 / fname
        return gff.exists() and gff.with_name(gff.name + ".tbi").exists()
    return (_CACHE / res.sha256 / "indexed.gff3.gz.tbi").exists()  # local build


def _all_grove_gg(name: str) -> Path:
    """Path to the lazily-built whole-genome `.gg` (may not exist yet)."""
    return _CACHE / "groves" / f"{RESOURCES[name].sha256}.{_GROVE_SCHEMA}" / "_all.gg"


def ensure_all_grove(name: str) -> Path:
    """Cache the whole-genome grove (`.gg`) if absent, returning its path.

    Prefers the **pinned prebuilt grove** (``grove_url``): a ~105 MB sha-verified download
    (seconds) instead of a local build. Falls back to building from the annotation
    (``build_grove(region="")`` → serialize) — minutes, and only when the resource declares no
    ``grove_layers``, because a local build reads only the annotation and cannot reproduce them.
    Either way it's cached; located queries never trigger this.

    Raises ``RuntimeError`` rather than returning a grove that is missing declared layers: the
    failure would otherwise be silent, since a query against an absent layer returns an empty
    result rather than an error.
    """
    gg = _all_grove_gg(name)
    if gg.exists():
        return gg
    res = RESOURCES[name]
    gg.parent.mkdir(parents=True, exist_ok=True)
    if res.grove_url:  # download the pinned .gg (fast, reproducible)
        # Short label: it repeats on every progress line, and a line wider than the terminal
        # wraps — after which `\r` cannot overwrite it cleanly. The contents are named once in
        # the header message instead (see `cli._prepare`).
        out = _download(res.grove_url, res.grove_sha256, gg, label="grove")
        _prune_superseded_groves(name)  # only after a verified fetch — see the docstring
        return out
    if res.grove_layers:
        raise RuntimeError(
            f"{name}: no `grove_url` is pinned, and a local build cannot reproduce this grove's "
            f"extra layer(s): {', '.join(res.grove_layers)}. Building from the annotation alone "
            "would return a structurally different grove whose queries come back empty instead of "
            "failing. Pin the prebuilt artifact, or clear `grove_layers` if it is genuinely "
            "annotation-only."
        )
    from genogrove_canopy.gff import build_grove  # local fallback — annotation-only resource

    tmp = gg.with_name(gg.name + ".tmp")
    build_grove(indexed_path(name), region="").serialize(str(tmp))
    tmp.replace(gg)
    _prune_superseded_groves(name)
    return gg


def _prune_superseded_groves(name: str) -> int:
    """Delete cached grove directories for ``name`` under a *different* ``_GROVE_SCHEMA``.

    Only ever called after a successful fetch or build, so what it removes is already superseded by
    a verified artifact. Without this every schema bump leaks a full grove: one real install carried
    988 MB across `.1` and `.2` before anyone noticed, meaning the bump *before* that had already
    leaked one silently.

    Matches on the resource's own sha256 prefix, so it cannot touch another resource's cache, and
    keeps both current directories (the grove and its `.shards` sibling). Returns how many it
    removed.
    """
    root = _CACHE / "groves"
    if not root.is_dir():
        return 0
    sha = RESOURCES[name].sha256
    keep = {f"{sha}.{_GROVE_SCHEMA}", f"{sha}.{_GROVE_SCHEMA}.shards"}
    removed = 0
    for d in root.glob(f"{sha}.*"):
        if d.name in keep or not d.is_dir():
            continue
        shutil.rmtree(d, ignore_errors=True)
        removed += 1
    return removed


def grove_view(name: str):
    """Open the whole-genome grove as a lazy ``GroveView`` (downloads the pinned `.gg`
    on first use, else builds it). Serves both located and genome-wide queries — pages in
    only the blocks a query touches; no whole-grove load, no per-query rebuild."""
    import pygenogrove as pg

    return pg.GroveView.open(str(ensure_all_grove(name)))


# Bump when the grove's content or payload model changes, so a stale `.gg` is re-fetched/rebuilt
# rather than silently served. v2 = a feature's `biotype` is its own (transcript_type on a
# transcript), exons carry none, and non-hierarchy types keep their column-9 attributes + `source`.
# v3 = the pinned artifact is the **unified** grove (backbone + ENCODE cCREs in one build).
# v4 = exons are external (graph-only) keys, not indexed in the B+ tree — reach one via its
# transcript's first_exon/next chain, never via intersect().
# v5 = transcripts are also external (reach via their gene's `contains` edge); only gene is
# indexed among the GENCODE hierarchy.
# v6 = v5 without the `intergenic_region` gap nodes it had briefly baked in for the SV layer;
# the pinned artifact carries no layer-specific structure (see `layers/sv.py`).
# v7 = `stop_codon_redefined_as_selenocysteine` (the selenocysteine-readthrough analogue of
# `stop_codon`) is now dropped like every other coding-structure annotation, instead of falling
# through to the generic foreign-layer path and getting indexed with its raw GENCODE attributes
# — 155 stray nodes genome-wide, found by checking what the indexed `type`s actually were.
# v8 = same content as v7, re-serialized under the pinned pygenogrove 0.9.0 (stream format 0.3).
# The v7 artifact was written by a stale 0.7.4 install (format 0.2), which 0.9.0 rejects with
# "bad magic" — and a cached copy is never re-fetched, so the only way to heal an install that
# already held it is a new schema directory.
#
# The bump is load-bearing, not cosmetic: this value is part of the cache directory
# (`_grove_dir`), but the rest of that key is the *annotation's* sha256 — which did not change when
# `grove_url` was re-pinned. Without the bump, anyone holding a cached `_all.gg` from the old
# GENCODE-only Zenodo artifact would keep being served it forever, since `ensure_all_grove` returns
# early on an existing path and never re-checks the URL.
_GROVE_SCHEMA = "8"


def _grove_dir(name: str) -> Path:
    """Directory for the **structure-only sharded index** — deliberately a *sibling* of the
    directory holding the pinned whole-genome grove, not the same one.

    They used to share both the directory and the `_all.gg` filename while being produced by
    different code from different sources: `ensure_all_grove` downloads the pinned artifact
    (annotation + its baked layers), while `grove_index` builds shards from the annotation alone,
    filtered to ``{gene, transcript, exon}``. Sharing meant `grove_index` treated a downloaded
    grove as proof its own shards existed, and `load_grove`'s rebuild-on-failure could delete the
    verified pinned artifact and silently replace it with a structure-only one.
    """
    return _CACHE / "groves" / f"{RESOURCES[name].sha256}.{_GROVE_SCHEMA}.shards"


def grove_index(name: str) -> tuple[dict[str, str], str]:
    """Resolve ``name`` to its sharded grove index, building + caching once.

    Returns ``({seqid: shard_path}, all_path)``: one serialized `.gg` per
    chromosome plus a whole-genome ``_all.gg``. A query deserializes only the
    shard(s) for the chromosome(s) it touches (fast, low-memory); ``_all`` is the
    whole-genome grove for genome-wide or cross-chromosome queries. Built in one
    streaming pass on first use (``genogrove_canopy.gff.write_sharded_groves``) and cached under
    ``<cache>/groves/<sha>.<schema>.shards/``; bump ``_GROVE_SCHEMA`` for model changes.

    **Structure-only, and not a substitute for the shipped grove.** The build filters to
    ``{gene, transcript, exon}``, so any layer the pinned artifact carries (the ENCODE cCREs)
    is absent here by construction. Nothing in the query path uses this — ``_grove_context``
    goes through ``ensure_all_grove``. Keep the two apart: see ``_grove_dir``.
    """
    from genogrove_canopy.gff import write_sharded_groves

    d = _grove_dir(name)
    if not (d / "_all.gg").exists():
        src = resolve(name)  # downloads + sha256-verifies the source once
        tmp = d.with_name(d.name + ".tmp")
        shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True, exist_ok=True)
        write_sharded_groves(src, tmp, types={"gene", "transcript", "exon"})
        shutil.rmtree(d, ignore_errors=True)
        tmp.replace(d)  # swap the finished index in atomically
    shards = {p.stem: str(p) for p in d.glob("*.gg") if p.name != "_all.gg"}
    return shards, str(d / "_all.gg")


def grove_path(name: str) -> Path:
    """The whole-genome `.gg` for ``name`` (the ``_all`` grove); builds the index if absent."""
    return Path(grove_index(name)[1])


def is_grove_cached(name: str) -> bool:
    """True if ``name``'s grove index is already built (so resolving it won't rebuild)."""
    return (_grove_dir(name) / "_all.gg").exists()


def load_grove(name: str):
    """Deserialize the **structure-only** sharded grove for ``name`` (cached). Self-heals once if
    the cached index won't deserialize.

    The rebuild is destructive, so it is scoped to the shard directory (``_grove_dir``) — never the
    directory holding the pinned artifact. It previously shared that directory, which meant a single
    transient deserialize failure could delete the sha-verified download and replace it with a
    structure-only rebuild, permanently and without a message.

    Like the rest of the shard path this is not used by the query path; prefer
    ``ensure_all_grove``/``grove_view``, which serve the pinned grove including its layers.
    """
    import pygenogrove as pg

    try:
        return pg.Grove.deserialize(str(grove_path(name)))
    except Exception:
        if _grove_dir(name) == _all_grove_gg(name).parent:  # belt-and-braces: never fire
            raise RuntimeError(
                f"refusing to rebuild {name}: the shard directory is the same as the pinned "
                "grove's, so the rebuild would delete a sha-verified artifact"
            ) from None
        shutil.rmtree(_grove_dir(name), ignore_errors=True)  # nuke the shard index -> rebuild
        return pg.Grove.deserialize(str(grove_path(name)))


# --------------------------------------------------------------------------- #
# ENCODE-rE2G enhancer→gene predictions — two lazy axes (see the design memo):
#   1. biosample: the catalog is metadata; only a requested biosample's edge BED is
#      ever fetched (1 of ~1,460), the rest cost a catalog row.
#   2. region: a fetched BED is bgzip+tabix-indexed once, then a query reads only the
#      rows overlapping its locus — never the whole ~90k-edge file.
# So a query's footprint = (biosamples it names) × (rows in the region it names).
# --------------------------------------------------------------------------- #


# The rE2G thresholded element-gene-links BED has ~56 columns (mostly model-internal
# `.Feature` inputs); we keep only these, selected BY HEADER NAME — the file's own
# `#`-prefixed header line — not by position, because `Score` is the LAST column and the
# tail count/order varies by rE2G version. `class` is a positional label (col 5); `Score`
# (col 56) is the model's calibrated enhancer→gene confidence. Everything else is either
# derivable from the grove (distance, class) or a model input we don't store.

# Catalog of all ENCODE-rE2G prediction annotations, generated by
# tools/fetch_re2g_catalog.py (one row per biosample; ships with the package).
_RE2G_CATALOG = Path(__file__).parent / "data" / "encode_re2g_catalog.tsv"


def re2g_catalog() -> list[dict[str, str]]:
    """The ENCODE-rE2G biosample catalog (metadata only — no data fetched).

    Each entry maps a biosample to its annotation ``accession``; the agent scopes which
    biosample(s) a query needs (axis 1). Raises ``FileNotFoundError`` with a pointer to
    the fetcher if the catalog hasn't been generated yet.
    """
    if not _RE2G_CATALOG.exists():
        raise FileNotFoundError(
            f"{_RE2G_CATALOG} missing — run `python tools/fetch_re2g_catalog.py` once "
            "to generate the ENCODE-rE2G catalog."
        )
    import csv

    with _RE2G_CATALOG.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# Ontology-id prefix -> the query axis it represents (UBERON = anatomy/tissue, CL = cell
# type, EFO = cell line / disease, NTR = ENCODE novel term). Lets the agent filter cohorts
# by axis without any external ontology file.
_ONTOLOGY_AXIS = {"UBERON": "tissue", "CL": "cell type", "CLO": "cell line",
                  "EFO": "cell line", "NTR": "novel term"}



# --------------------------------------------------------------------------- #
# ENCODE-rE2G index bundle — the enhancer→gene layer, pinned per file.
#
# 369 cohorts x (byEnhancer, byTargetGene) x (.tsv.gz, .tbi) = 1,476 files, ~3 GB. Too many for
# `RESOURCES`, which is a hand-written catalog of named datasets, so the checksums live in a
# manifest shipped with the package and the URL is derived from one pinned commit.
#
# These are the *derived* bgzip+tabix files a query reads, not the raw ENCODE BEDs: pinning them
# means a user needs no htslib and no local sort/index step, and gets exactly the bytes the
# published results were computed from.
# --------------------------------------------------------------------------- #

#: Immutable commit holding the index bundle. Same rule as every other pin — never a branch.
RE2G_INDEX_COMMIT = "38283ced34d85edfe6ee2b936ca7f17db535bd84"
RE2G_INDEX_MANIFEST = Path(__file__).resolve().parent / "data" / "re2g_index.manifest.tsv"


@lru_cache(maxsize=1)
def re2g_index_manifest() -> dict[str, str]:
    """``{filename: sha256}`` for every file in the pinned rE2G index bundle."""
    out = {}
    with open(RE2G_INDEX_MANIFEST, encoding="utf-8") as fh:
        next(fh)  # header
        for line in fh:
            name, digest, _size = line.rstrip("\n").split("\t")
            out[name] = digest
    return out


def re2g_index_file(filename: str) -> Path:
    """Resolve one file of the rE2G index bundle, downloading + verifying it on first use.

    Per file rather than per bundle: a question about one cohort should not pull 3 GB. The four
    files a cohort needs are ~8 MB together.
    """
    dest = _CACHE / "re2g_index" / filename
    if dest.exists():
        return dest
    manifest = re2g_index_manifest()
    if filename not in manifest:
        raise KeyError(f"{filename!r} is not in the pinned rE2G index manifest")
    if not RE2G_INDEX_COMMIT:
        raise RuntimeError(
            "no rE2G index commit is pinned, so the enhancer layer cannot be fetched. "
            "Set RE2G_INDEX_COMMIT to the genogrove/canopy commit holding `re2g/`."
        )
    url = (f"https://huggingface.co/datasets/genogrove/canopy/resolve/"
           f"{RE2G_INDEX_COMMIT}/re2g/{filename}")
    return _download(url, manifest[filename], dest, label=f"rE2G {filename}")


def re2g_cohorts() -> list[dict]:
    """The catalog grouped into **cohorts** — one per biosample (its ontology id), folding
    the replicate accessions together. This is the request→biosample association layer: an
    agent picks a cohort from this list by ``name`` / ``axis`` / ``type`` (grounded — only
    declared cohorts exist), and the cohort's ``accessions`` then drive a lazy per-cohort
    grove build (merge replicates → edge support ``n``). Derived from the catalog, so it
    ships nothing new.

    Each cohort: ``ontology_id``, ``name`` (biosample term), ``axis`` (from the ontology
    prefix), ``type`` (biosample_type), ``n_replicates``, ``accessions``. Sorted by replicate
    count (most-replicated first), then name.
    """
    groups: dict[str, dict] = {}
    for e in re2g_catalog():
        g = groups.setdefault(e["biosample_id"], {
            "ontology_id": e["biosample_id"], "name": e["biosample_term"],
            "axis": _ONTOLOGY_AXIS.get(e["biosample_id"].split(":")[0], "other"),
            "type": e["biosample_type"], "accessions": [],
        })
        g["accessions"].append(e["accession"])
    for g in groups.values():
        g["n_replicates"] = len(g["accessions"])
    return sorted(groups.values(), key=lambda g: (-g["n_replicates"], g["name"]))


