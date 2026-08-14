# SPDX-License-Identifier: GPL-3.0-or-later
"""Guard the Level 2 reproducibility invariant: the pinned ``pygenogrove`` build
must be described identically in four places — the ``==`` dependency pin and the
``[tool.uv.sources]`` ``rev`` in ``pyproject.toml``, ``resources.PYGENOGROVE``, and
the target version stated in the ``prompts/system.md`` header.

If any of the four drifts, a run records a build it was not actually made
against — or, for the prompt, generates code against a surface the pinned build
may not have. This test mechanizes the manual "all four agree" QC check so drift
fails CI instead of relying on review.

``pyproject.toml`` is parsed with regexes rather than ``tomllib`` so the test
runs on the supported floor (py3.9, where ``tomllib`` does not exist) in the
ephemeral ``--with pytest`` CI env (no ``tomli`` available either).
"""

from __future__ import annotations

import re
from pathlib import Path

from genogrove_canopy import resources
from genogrove_canopy.resources import PYGENOGROVE

PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"
# Derived from the package rather than the repo layout, so it follows a package rename (the
# import name changed in #6) and points at the file `llm.py` actually loads the prompt from.
SYSTEM_MD = Path(resources.__file__).resolve().parent / "prompts" / "system.md"


def _pyproject_text() -> str:
    return PYPROJECT.read_text(encoding="utf-8")


def test_dependency_pin_matches_registry() -> None:
    """The ``pygenogrove==<version>`` dependency pin matches ``PYGENOGROVE.version``."""
    text = _pyproject_text()
    m = re.search(r'pygenogrove==(?P<version>[0-9][^"\'\s]*)', text)
    assert m, "no `pygenogrove==<version>` dependency pin found in pyproject.toml"
    assert m.group("version") == PYGENOGROVE.version


def test_uv_source_rev_matches_registry() -> None:
    """The ``[tool.uv.sources]`` ``rev`` matches ``PYGENOGROVE.git_rev`` exactly."""
    text = _pyproject_text()
    m = re.search(
        r'pygenogrove\s*=\s*\{[^}]*\brev\s*=\s*"(?P<rev>[0-9a-f]{40})"',
        text,
    )
    assert m, "no `pygenogrove = { ... rev = \"<40-hex>\" }` source found in pyproject.toml"
    assert m.group("rev") == PYGENOGROVE.git_rev


def test_system_prompt_target_matches_registry() -> None:
    """The version ``prompts/system.md`` names as its target matches ``PYGENOGROVE.version``.

    ``system.md`` *is* the system prompt, so this header ships to the model on every
    request — a stale version tells it the documented API surface targets a build that
    isn't the one its code will run against. Nothing enforced this before, and the
    header silently fell three releases behind (canopy#7).
    """
    text = SYSTEM_MD.read_text(encoding="utf-8")
    # trailing `[0-9A-Za-z]` so the sentence's final "." isn't captured as part of the version
    m = re.search(
        r"Current target:\s*pygenogrove\s+(?P<version>[0-9](?:[0-9A-Za-z.\-+]*[0-9A-Za-z])?)",
        text,
    )
    assert m, "no `Current target: pygenogrove <version>` line found in prompts/system.md"
    assert m.group("version") == PYGENOGROVE.version, (
        f"prompts/system.md targets pygenogrove {m.group('version')} but the pinned build is "
        f"{PYGENOGROVE.version} — update the header in the same commit as the pin"
    )


def test_every_pinned_artifact_is_immutable() -> None:
    """**Every** pinned artifact — annotation, index, grove — must be checksummed and immutable.

    The build pin has had a drift guard since #2; the *data* pins never did, which is how a grove
    built under an older payload model stayed pinned across three releases (#9). A movable ref is
    the failure mode that matters: `resolve/main` on Hugging Face or a `/tree/` URL keeps resolving
    after the bytes behind it change, so the checksum starts failing on a fetch that used to work —
    or worse, someone updates the sha to match and the pin silently means something else.

    This covered only `grove_url` until #14 gave `encode.ccre.v4` real Hugging Face URLs, which the
    guard then did not look at: it skipped every resource without a grove. Now it walks all three
    (url, index_url, grove_url) with their matching checksums.
    """
    movable = ("/resolve/main/", "/resolve/master/", "/tree/", "/blob/main/", "/raw/main/")
    for name, res in resources.RESOURCES.items():
        for field, url, digest in (
            ("url", res.url, res.sha256),
            ("index_url", res.index_url, res.index_sha256),
            ("grove_url", res.grove_url, res.grove_sha256),
        ):
            if not url:
                continue
            assert re.fullmatch(r"[0-9a-f]{64}", digest), (
                f"{name}.{field}: needs a full 64-hex sha256, got {digest!r}"
            )
            for ref in movable:
                assert ref not in url, (
                    f"{name}.{field} pins the movable ref {ref!r} — use an immutable commit sha so "
                    f"the pinned bytes cannot change underneath the checksum ({url})"
                )
            # A Hugging Face `resolve/<ref>/…` must name a 40-hex commit, not a branch or tag.
            m = re.search(r"huggingface\.co/datasets/[^/]+/[^/]+/resolve/(?P<ref>[^/]+)/", url)
            if m:
                assert re.fullmatch(r"[0-9a-f]{40}", m.group("ref")), (
                    f"{name}.{field} resolves ref {m.group('ref')!r}, which is not an immutable "
                    "commit sha"
                )


def test_no_resource_url_is_a_placeholder() -> None:
    """**Every** catalog URL must be reachable — no placeholders anywhere.

    Until the cCRE pair was uploaded, `encode.ccre.v4` carried `pending-upload.invalid` URLs and
    those 32 MB existed only in one developer's cache: a fresh install could not bake the grove
    (#9), and even after the layer moved into the grove, rebuilding it depended on that machine.
    The exemption this guard used to carry is gone, so any new placeholder fails here.
    """
    for name, res in resources.RESOURCES.items():
        urls = [u for u in (res.url, res.index_url, res.grove_url) if u]
        placeholders = [u for u in urls if ".invalid" in u or "pending-upload" in u]
        assert not placeholders, (
            f"{name}: placeholder URL(s) {placeholders} would fail on any machine without a "
            "pre-seeded cache"
        )


def test_pin_is_an_immutable_commit() -> None:
    """The pin must be a full 40-char commit SHA, not a movable branch/tag."""
    assert re.fullmatch(r"[0-9a-f]{40}", PYGENOGROVE.git_rev), (
        "PYGENOGROVE.git_rev must be a full immutable commit SHA "
        f"(got {PYGENOGROVE.git_rev!r})"
    )


def test_declared_grove_layers_cannot_be_silently_dropped(monkeypatch, tmp_path) -> None:
    """A resource declaring `grove_layers` must refuse a local build, not return a lesser grove.

    A local build reads only the annotation, so it cannot reproduce a baked layer. Returning one
    anyway is the worst failure shape available: queries against the missing layer come back
    **empty rather than failing**, so nothing upstream notices. Before #11 this could not happen
    because the bake added the layer to whatever grove came back; removing the bake removed that
    safety net, so the refusal replaces it.
    """
    import dataclasses

    res = resources.RESOURCES["gencode.human"]
    assert res.grove_layers, "gencode.human must declare the layers its pinned grove carries"

    # Point the cache at an empty dir: `ensure_all_grove` returns early when the grove is already
    # cached, so a warm cache never reaches the branch under test.
    monkeypatch.setattr(resources, "_CACHE", tmp_path / "cache")
    monkeypatch.setitem(
        resources.RESOURCES, "gencode.human",
        dataclasses.replace(res, grove_url="", grove_sha256=""),
    )
    try:
        resources.ensure_all_grove("gencode.human")
    except RuntimeError as exc:
        assert "layer" in str(exc), f"refused, but not for the layer reason: {exc}"
    else:
        raise AssertionError(
            "ensure_all_grove built a grove locally despite declared grove_layers — queries "
            "against the absent layer would silently return empty"
        )


def test_shard_index_is_not_the_pinned_grove_directory() -> None:
    """The structure-only shard index must live in a *different* directory from the pinned grove.

    They shared a directory and the `_all.gg` filename while being built by different code from
    different sources. Two consequences, both silent: `grove_index` treated a downloaded grove as
    proof its own shards existed, and `load_grove`'s rebuild-on-failure deleted the whole directory
    — destroying the sha-verified artifact and replacing it with a `{gene,transcript,exon}`-filtered
    rebuild that has no cCREs.
    """
    for name in resources.RESOURCES:
        shard_dir = resources._grove_dir(name)
        pinned_dir = resources._all_grove_gg(name).parent
        assert shard_dir != pinned_dir, (
            f"{name}: shard index and pinned grove share {shard_dir} — a shard rebuild would "
            "delete the pinned artifact"
        )
        assert pinned_dir not in shard_dir.parents, (
            f"{name}: shard index nests inside the pinned grove's directory ({shard_dir})"
        )


def test_prune_removes_only_superseded_grove_dirs(monkeypatch, tmp_path) -> None:
    """Pruning deletes stale-schema dirs for this resource and touches nothing else.

    Every `_GROVE_SCHEMA` bump used to leak a whole grove: one real install carried 988 MB across
    `.1` and `.2` before anyone looked, which means the bump before that had already leaked one.
    """
    monkeypatch.setattr(resources, "_CACHE", tmp_path)
    sha = resources.RESOURCES["gencode.human"].sha256
    cur = resources._GROVE_SCHEMA
    root = tmp_path / "groves"
    names = [f"{sha}.{cur}", f"{sha}.{cur}.shards",          # current — keep
             f"{sha}.1", f"{sha}.2", f"{sha}.1.shards",      # superseded — remove
             "0" * 64 + f".{cur}"]                            # another resource — never touch
    for n in names:
        (root / n).mkdir(parents=True)

    removed = resources._prune_superseded_groves("gencode.human")
    left = sorted(p.name for p in root.iterdir())
    assert removed == 3
    assert left == sorted([f"{sha}.{cur}", f"{sha}.{cur}.shards", "0" * 64 + f".{cur}"])


def test_cache_location_is_overridable(monkeypatch) -> None:
    """`GENOGROVE_CANOPY_CACHE` redirects the cache, so a cold run needs no deletion.

    Read at import, so this reloads the module rather than setting the variable and hoping.
    """
    import importlib

    monkeypatch.setenv("GENOGROVE_CANOPY_CACHE", "/tmp/canopy-cache-probe")
    reloaded = importlib.reload(resources)
    try:
        assert str(reloaded._CACHE) == "/tmp/canopy-cache-probe"
    finally:
        monkeypatch.delenv("GENOGROVE_CANOPY_CACHE")
        importlib.reload(resources)


def test_re2g_digest_is_recorded_and_readable(monkeypatch, tmp_path) -> None:
    """rE2G is the one layer fetched without a pinned checksum, so record what arrived.

    Pinning all 369 cohorts is a separate question (#20); this makes drift *detectable* after
    the fact rather than invisible, which is the part that costs nothing.
    """
    monkeypatch.setattr(resources, "_CACHE", tmp_path)
    raw = tmp_path / "raw.bed.gz"
    raw.write_bytes(b"chr1\t1\t2\tgene\n")

    digest = resources._record_re2g_digest("ENCSR000AAA", raw)
    assert len(digest) == 64
    assert resources.re2g_digest("ENCSR000AAA") == digest


def test_re2g_provenance_is_honest_about_what_it_does_not_know(monkeypatch, tmp_path) -> None:
    """A cohort cached before digests were recorded reports an empty string, not a guess.

    The raw BED is deleted after indexing, so there is nothing to back-fill from — the indexed
    file has been re-sorted and re-compressed and would hash differently.
    """
    monkeypatch.setattr(resources, "_CACHE", tmp_path)
    raw = tmp_path / "raw.bed.gz"
    raw.write_bytes(b"x")
    known = resources._record_re2g_digest("ENCSR000KNW", raw)

    prov = resources.re2g_provenance(["ENCSR000KNW", "ENCSR000UNK"])
    assert prov["ENCSR000KNW"] == known
    assert prov["ENCSR000UNK"] == "", "an unknown digest must not be reported as anything else"
