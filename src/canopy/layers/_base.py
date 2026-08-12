# SPDX-License-Identifier: GPL-3.0-or-later
"""The Layer contract — one descriptor per omic layer.

A ``Layer`` carries the metadata the **agent** reads to decide whether a question implicates
this layer and what its query code will see (``axis``/``kind``/``when``/``schema``), plus the
uniform host-side ``attach`` op that materialises fetched records as grove nodes/edges. Each
layer module (``canopy.layers.ccres``, ``canopy.layers.enhancers``, …) defines its own
``fetch`` shaped to its need and exports a module-level ``LAYER = Layer(...)``; the registry in
``canopy.layers`` collects them. See ``docs/showcase-design.md`` §4 (taxonomy) and §6
(on-demand resolution).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass


@dataclass(frozen=True)
class Layer:
    """Describes one loadable layer to the agent and the host.

    ``axis``   — coordinate axis it lives on: ``genomic`` | ``transcript`` | ``protein`` | ``cross``.
    ``kind``   — taxonomy: ``node`` (new intervals) | ``edge`` (relations) | ``payload`` (values).
    ``when``   — one line: when a question implicates this layer (drives §6 selection).
    ``schema`` — the node/edge ``type`` + payload keys the generated query code will see.
    ``attach`` — ``(grove, records) -> int``: insert the fetched records; returns the count.
    """

    name: str
    axis: str
    kind: str
    title: str
    when: str
    schema: str
    attach: Callable[..., int]
