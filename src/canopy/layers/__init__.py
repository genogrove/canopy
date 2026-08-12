# SPDX-License-Identifier: GPL-3.0-or-later
"""Layer registry — the catalogue the agent selects from.

Each layer module exports a ``LAYER`` descriptor (see ``_base.Layer``); this package collects
them into ``REGISTRY`` and renders the agent-facing catalogue. Adding a layer = adding a module
here and listing it below — the system prompt updates from the registry, so there is one source
of truth for "what data exists and how the query code sees it" (``docs/showcase-design.md`` §6).
"""

from __future__ import annotations

from canopy.layers import ccres, enhancers
from canopy.layers._base import Layer

# One entry per layer module. Order is the order shown to the agent.
REGISTRY: dict[str, Layer] = {m.LAYER.name: m.LAYER for m in (ccres, enhancers)}


def catalogue_block(names: list[str] | None = None) -> str:
    """The agent-facing layer catalogue as prompt text. ``names`` limits it to the implicated
    layers (§6 on-demand injection); ``None`` lists every registered layer."""
    layers = REGISTRY.values() if names is None else (REGISTRY[n] for n in names)
    return "\n".join(
        f"- **{l.name}** ({l.axis} axis · {l.kind}) — {l.title}. "
        f"Load when {l.when}. Query code sees: {l.schema}"
        for l in layers
    )


__all__ = ["REGISTRY", "Layer", "catalogue_block", "ccres", "enhancers"]
