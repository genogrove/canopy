# SPDX-License-Identifier: GPL-3.0-or-later
"""Guard the Level 2 reproducibility invariant: the pinned ``pygenogrove`` build
must be described identically in three places — the ``==`` dependency pin and the
``[tool.uv.sources]`` ``rev`` in ``pyproject.toml``, and ``resources.PYGENOGROVE``.

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

PYPROJECT = Path(__file__).resolve().parent.parent / "pyproject.toml"


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


def test_pin_is_an_immutable_commit() -> None:
    """The pin must be a full 40-char commit SHA, not a movable branch/tag."""
    assert re.fullmatch(r"[0-9a-f]{40}", PYGENOGROVE.git_rev), (
        "PYGENOGROVE.git_rev must be a full immutable commit SHA "
        f"(got {PYGENOGROVE.git_rev!r})"
    )
