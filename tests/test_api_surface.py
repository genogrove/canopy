# SPDX-License-Identifier: GPL-3.0-or-later
"""Guard the codegen contract against the *pinned build's actual API surface*.

``tests/test_resources_pins.py`` guards the version **string** the prompt claims to target. That
is a pure text comparison and catches a stale label. It cannot catch the failure that actually
breaks a run: the prompt advertising a method the pinned ``pygenogrove`` does not have. Generated
code then raises ``AttributeError`` inside the sandbox and the user sees a broken answer rather
than a build error — no checksum, pin, or string match notices.

It has gone wrong in both directions already:

* **Overstating** — ``457feaa`` removed ``get_edge_list`` from the prompt because it was *not* on
  ``GroveView``. So the prompt had been advertising a method the build lacked.
* **Understating** — ``get_edge_list`` *is* on ``GroveView`` in 0.7.4, so that removal is now
  itself stale and the prompt hides available capability. Harmless, but it shows the drift is
  bidirectional and was only ever caught by hand.

These tests need the real bindings, so they ``importorskip`` and are **skipped by the default CI
job**, which runs ``uv run --no-project`` and never installs ``pygenogrove``. The path-filtered
``api-surface`` workflow does a real ``uv sync`` and is what makes them binding — see
``.github/workflows/api-surface.yml`` and canopy#10.
"""

from __future__ import annotations

import pytest

from genogrove_canopy.cli import QUERY_SURFACE

pg = pytest.importorskip("pygenogrove", reason="needs the built bindings (see api-surface workflow)")


def test_advertised_query_surface_exists_on_groveview() -> None:
    """Every method the prompt advertises must exist on the pinned build's ``GroveView``.

    This is the assertion that would have caught the pre-457feaa state, where the prompt named a
    method the build did not expose.
    """
    real = {n for n in dir(pg.GroveView) if not n.startswith("_")}
    missing = [m for m in QUERY_SURFACE if m not in real]
    assert not missing, (
        f"prompts advertise {missing} on GroveView, which pygenogrove does not expose. Generated "
        f"code calling them raises AttributeError in the sandbox. Available: {sorted(real)}"
    )


def test_advertised_surface_is_actually_read_only() -> None:
    """The prompt calls this surface "Query-only", so none of it may mutate the grove.

    ``GroveView`` is the read-side handle, but the claim is worth pinning: if a future build moved
    a mutating method onto it, the prompt would be licensing generated code to write to a
    sha-verified artifact the host owns.
    """
    mutating = ("insert", "insert_bulk", "insert_sorted", "add_edge", "remove_edge", "remove_key",
                "serialize", "compact", "clear_graph")
    overlap = [m for m in QUERY_SURFACE if m in mutating]
    assert not overlap, f"QUERY_SURFACE advertises mutating method(s) {overlap} as query-only"


def test_query_surface_is_not_empty_and_is_used_by_the_prompt() -> None:
    """`QUERY_SURFACE` must be non-empty and must actually reach the model.

    Guards the orphan case: a constant that nothing interpolates would let this whole file pass
    while the prompt says something else entirely.
    """
    assert QUERY_SURFACE, "QUERY_SURFACE is empty — the prompt would advertise no query surface"

    from pathlib import Path

    from genogrove_canopy import cli

    src = Path(cli.__file__).read_text(encoding="utf-8")
    assert "QUERY_SURFACE" in src.split("QUERY_SURFACE = ", 1)[1], (
        "QUERY_SURFACE is defined but never referenced — the resources block must interpolate it, "
        "or this guard is checking a constant the model never sees"
    )
