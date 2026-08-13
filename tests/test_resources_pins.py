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


def test_pinned_grove_artifact_is_immutable() -> None:
    """Every pinned grove artifact must be a checksummed, immutable reference.

    The build pin has had a drift guard since #2; the *data* pins never did, which is how a
    grove artifact built under an older payload model stayed pinned across three releases
    (#9). A movable ref is the failure mode that matters: `resolve/main` on Hugging Face or a
    `/tree/` URL keeps resolving after the file behind it changes, so the sha256 check starts
    failing on a fetch that used to work — or worse, the sha is updated to match and the pin
    silently now means something else.
    """
    movable = ("/resolve/main/", "/resolve/master/", "/tree/", "/blob/main/", "/raw/main/")
    for name, res in resources.RESOURCES.items():
        if not res.grove_url:
            continue
        assert re.fullmatch(r"[0-9a-f]{64}", res.grove_sha256), (
            f"{name}: grove_sha256 must be a full 64-hex sha256 (got {res.grove_sha256!r})"
        )
        for ref in movable:
            assert ref not in res.grove_url, (
                f"{name}: grove_url pins the movable ref {ref!r} — use an immutable commit sha "
                f"so the pinned bytes cannot change underneath the checksum ({res.grove_url})"
            )
        # A Hugging Face `resolve/<ref>/…` must name a 40-hex commit, not a branch or tag.
        m = re.search(r"huggingface\.co/datasets/[^/]+/[^/]+/resolve/(?P<ref>[^/]+)/", res.grove_url)
        if m:
            assert re.fullmatch(r"[0-9a-f]{40}", m.group("ref")), (
                f"{name}: grove_url resolves ref {m.group('ref')!r}, which is not an immutable "
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
