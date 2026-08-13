# SPDX-License-Identifier: GPL-3.0-or-later
"""Guard the codegen contract against the pinned build's **actual** API surface.

``test_resources_pins.py`` compares version *strings*, which catches a stale label but not the
failure that breaks a run: the prompt advertising a method the pinned ``pygenogrove`` does not
have. Generated code then raises ``AttributeError`` inside the sandbox and the user sees a broken
answer rather than an error anyone can act on.

It has gone wrong both ways. ``457feaa`` removed ``get_edge_list`` from the prompt because it was
*not* on ``GroveView`` — so the prompt had been advertising a method the build lacked. And
``get_edge_list`` *is* on ``GroveView`` in 0.7.4, so that removal is now stale in the other
direction. Both were caught by hand.

This module needs the real bindings and so is skipped wherever they are absent — including the
default CI job, which runs ``uv run --no-project``. The path-filtered ``api-surface`` workflow does
a real ``uv sync`` and is what makes it binding. Checks that do *not* need the bindings live in
``test_codegen_contract.py`` so they run everywhere; keep it that way.
"""

from __future__ import annotations

import pytest

from genogrove_canopy.cli import QUERY_SURFACE

pg = pytest.importorskip("pygenogrove", reason="needs the built bindings (see api-surface workflow)")


def test_advertised_query_surface_exists_on_groveview() -> None:
    """Every method the prompt advertises must exist on the pinned build's ``GroveView``.

    The assertion that would have caught the pre-457feaa state, where the prompt named a method
    the build did not expose.
    """
    real = {n for n in dir(pg.GroveView) if not n.startswith("_")}
    missing = [m for m in QUERY_SURFACE if m not in real]
    assert not missing, (
        f"prompts advertise {missing} on GroveView, which pygenogrove does not expose. Generated "
        f"code calling them raises AttributeError in the sandbox. Available: {sorted(real)}"
    )
