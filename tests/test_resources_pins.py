# SPDX-License-Identifier: GPL-3.0-or-later
"""Guard the Level 2 reproducibility invariant: the pinned ``pygenogrove`` build
must be described identically in four places — the ``==`` dependency pin and the
``[tool.uv.sources]`` ``rev`` in ``pyproject.toml``, ``resources.PYGENOGROVE``, and
the target version stated in the ``prompts/system.md`` header.

If any of the three drifts, a run records a build it was not actually made
against. This test mechanizes the manual "all three agree" QC check so drift
fails CI instead of relying on review.

``pyproject.toml`` is parsed with regexes rather than ``tomllib`` so the test
runs on the supported floor (py3.9, where ``tomllib`` does not exist) in the
ephemeral ``--with pytest`` CI env (no ``tomli`` available either).
"""

from __future__ import annotations

import re
from pathlib import Path

from genogrove_canopy.resources import PYGENOGROVE

_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = _ROOT / "pyproject.toml"
SYSTEM_MD = _ROOT / "src" / "genogrove_canopy" / "prompts" / "system.md"


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


def test_pin_is_an_immutable_commit() -> None:
    """The pin must be a full 40-char commit SHA, not a movable branch/tag."""
    assert re.fullmatch(r"[0-9a-f]{40}", PYGENOGROVE.git_rev), (
        "PYGENOGROVE.git_rev must be a full immutable commit SHA "
        f"(got {PYGENOGROVE.git_rev!r})"
    )
