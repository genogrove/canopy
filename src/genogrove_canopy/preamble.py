# SPDX-License-Identifier: GPL-3.0-or-later
"""The sandbox preamble: binds ``GROVE`` to an open, ready-to-query grove with the question's
layers attached, and memoises that grove across the queries of a warm worker.

Per-question layers (rE2G enhancer cohorts, PCAWG SV cohorts) are attached **inside the
sandbox**, because that is the only place the generated code can reach them: an in-memory
grove cannot cross the process boundary, and the records cannot be injected as a literal
(~34 MB of program text for one enhancer cohort alone). The layer's attach function is shipped
as source text (``inspect.getsource``), so the exact code the tests exercise is what runs.

The grove is memoised in ``_CANOPY_STATE`` (see ``sandbox.Worker``) **and grows**: a cohort
asked about later is attached onto the same grove, never rebuilt — so after a breast-cancer
question and then a prostate one the grove holds both, and "in prostate but not breast" is a
set operation on the ``byCohort`` map (enhancers) or the ``cohort`` field (SV edges). The
generated code is told which cohorts *this* question is about through ``COHORTS`` (rE2G ids)
and ``SV_COHORTS`` (PCAWG project codes) and must filter on them; the prompt says so.
ponytail: the grove only grows within a session — a few hundred MB per enhancer cohort, tens
of MB per SV cohort, no eviction. Add an LRU over the cohort sets if sessions that wander
across many tissues turn up.
"""

from __future__ import annotations

import inspect
import json

from genogrove_canopy.layers import enhancers, sv


def build(gg: str, cohort_links: dict[str, str] | None = None,
          sv_files: dict[str, str] | None = None) -> str:
    """The preamble for one question. ``cohort_links``: rE2G cohort id -> its links table
    (``enhancers.links_file``); ``sv_files``: PCAWG project code -> its cohort table
    (``sv.cohort_file``). With neither, ``GROVE`` is a lazy ``GroveView`` (~200 ms); with
    either, a mutable ``Grove`` from the memo, with the missing cohorts attached in turn."""
    gg_lit = json.dumps(gg)
    cohorts = sorted(cohort_links or {})
    svs = sorted(sv_files or {})
    if not cohorts and not svs:  # the worker hands _CANOPY_STATE to every query: hide it here too
        return (f"import pygenogrove as pg\nGROVE = pg.GroveView.open({gg_lit})\n"
                "COHORTS = []\nSV_COHORTS = []\nglobals().pop('_CANOPY_STATE', None)\n")

    return (
        "import pygenogrove as pg\n"
        f"{inspect.getsource(enhancers.attach_links)}\n"
        f"{inspect.getsource(sv.read_table)}\n"
        f"{inspect.getsource(sv.attach_tracked)}\n"
        "_state = globals().get('_CANOPY_STATE')\n"       # absent in one-shot `sandbox.run`
        "if _state is None:\n"
        "    _state = {}\n"
        f"if _state.get('gg') != {gg_lit}:\n"
        "    _state.clear()\n"   # drop a stale grove BEFORE deserializing the next one: both
        # live at once is ~1.8 GB + a ~2 GB deserialize peak against the 4 GiB sandbox cap
        f"    _state.update(gg={gg_lit}, grove=pg.Grove.deserialize({gg_lit}),\n"
        "                  nodes={}, cohorts=set(), sv_cohorts=set())\n"
        "_grove = _state['grove']\n"
        # A cohort is marked attached only after its attach returned, so an exception halfway
        # would leave a half-mutated grove that the next query attaches onto again: drop the
        # memo entirely on failure and let the next query rebuild.
        "try:\n"
        f"    for _c, _path in {json.dumps(sorted((cohort_links or {}).items()))}:\n"
        "        if _c not in _state['cohorts']:\n"        # attach onto the warm grove, once
        "            attach_links(_grove, _path, _c, _state['nodes'])\n"
        "            _state['cohorts'].add(_c)\n"
        f"    for _c, _path in {json.dumps(sorted((sv_files or {}).items()))}:\n"
        "        if _c not in _state['sv_cohorts']:\n"
        "            with open(_path) as _fh:\n"
        "                attach_tracked(_grove, read_table(_fh))\n"
        "            _state['sv_cohorts'].add(_c)\n"
        "except BaseException:\n"
        "    _state.clear()\n"
        "    raise\n"
        f"COHORTS = {json.dumps(cohorts)}\n"               # this question's cohorts, by layer key
        f"SV_COHORTS = {json.dumps(svs)}\n"
        # The memoised grove is a mutable `Grove` shared by every query of the session. Hand the
        # generated code a view that forwards only what a read-only `GroveView` has (plus the
        # count/lookup helpers), so a query cannot insert into it and quietly change the next
        # answer. Closures, not bound methods, so no `.__self__` leads back to the grove.
        "def _readonly(g):\n"
        "    _ok = {n for n in dir(pg.GroveView) if not n.startswith('_')} | {\n"
        "        'size', 'edge_count', 'vertex_count', 'external_vertex_count',\n"
        "        'indexed_vertex_count', 'vertex_count_with_edges', 'key_storage_size',\n"
        "        'has_edge', 'graph_empty'}\n"
        "    class _View:\n"
        "        def __getattr__(self, n):\n"
        "            if n in _ok:\n"
        "                return lambda *a, **k: getattr(g, n)(*a, **k)\n"
        "            raise AttributeError(f'GROVE is read-only: {n!r} is not a query method')\n"
        "        def __len__(self):\n"
        "            return len(g)\n"
        "    return _View()\n"
        "GROVE = _readonly(_grove)\n"
        # Drop the host-only names from the namespace the generated code runs in, so a query
        # won't by accident evict or swap the grove the next query in a warm session is given.
        # (`sys.modules['__main__']` can still reach it — the sandbox's documented residual risk.)
        "for _n in ('_CANOPY_STATE', '_state', '_grove', '_c', '_path', '_fh', '_readonly', '_n'):\n"
        "    globals().pop(_n, None)\n"
    )
