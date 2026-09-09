# SPDX-License-Identifier: GPL-3.0-or-later
"""The sandbox preamble: lazy view without cohorts; with cohorts a memoised mutable grove that
grows one cohort (enhancer or SV) at a time, a read-only view over it, and no host names left
in the namespace. Executed against stub bindings; the real-bindings runs are in
test_sandbox_pygenogrove.py."""
import pytest

from genogrove_canopy import preamble


def test_preamble_no_cohorts_opens_lazily():
    pre = preamble.build("/tmp/x.gg")
    assert pre == ('import pygenogrove as pg\n'
                    'GROVE = pg.GroveView.open("/tmp/x.gg")\nCOHORTS = []\nSV_COHORTS = []\n'
                    "globals().pop('_CANOPY_STATE', None)\n")
    compile(pre, "<preamble>", "exec")


def test_preamble_with_cohorts_attaches_each_onto_one_grove():
    pre = preamble.build(
        "/tmp/x.gg", {"EFO:0005726": "/tmp/a.tsv", "EFO:0009318": "/tmp/b.tsv"})
    assert 'grove=pg.Grove.deserialize("/tmp/x.gg")' in pre
    assert '[["EFO:0005726", "/tmp/a.tsv"], ["EFO:0009318", "/tmp/b.tsv"]]' in pre
    assert "attach_links(_grove, _path, _c, _state['nodes'])" in pre
    assert 'COHORTS = ["EFO:0005726", "EFO:0009318"]' in pre
    assert "GROVE = _readonly(_grove)" in pre
    assert "_CANOPY_STATE" in pre
    compile(pre, "<preamble>", "exec")


def _run_preamble(cohort_links, state, events, sv_files=None):
    """Execute the cohort preamble against stub bindings; returns the namespace generated code
    would see. `events` records deserialize / attach / clear calls in order."""

    class _Grove:
        n = 1
        def size(self):
            return self.n
        def __len__(self):
            return self.n

    class _GroveView:  # what the sandbox's read-only handle exposes
        def size(self): ...

    pg = type("pg", (), {
        "Grove": type("G", (), {"deserialize": staticmethod(
            lambda p: events.append("deser") or _Grove())}),
        "GroveView": _GroveView})
    pre = preamble.build("/tmp/x.gg", cohort_links, sv_files)
    body = pre.split("import pygenogrove as pg\n", 1)[1].replace(
        "        attach_links(_grove, _path, _c, _state['nodes'])",
        "        events.append(('attach', _c, id(_state['nodes'])))").replace(
        "        with open(_path) as _fh:\n            attach_tracked(_grove, read_table(_fh))",
        "        events.append(('sv', _c))")
    g = {"_CANOPY_STATE": state, "__builtins__": __builtins__, "pg": pg, "events": events}
    exec(body, g)
    return g


def test_preamble_grows_the_warm_grove_one_cohort_at_a_time():
    """Cohorts asked about later attach onto the memoised grove — never a rebuild — sharing one
    node cache so `attach_links` merges them; a question whose cohorts are all present attaches
    nothing; a different grove path drops the old state before deserializing the new grove."""

    class _State(dict):
        def clear(self):
            events.append("clear")
            super().clear()

    events, state = [], _State()
    _run_preamble({"C": "/tmp/c.tsv"}, state, events)
    assert events == ["clear", "deser", ("attach", "C", id(state["nodes"]))]
    assert state["cohorts"] == {"C"}

    events.clear()
    g = _run_preamble({"D": "/tmp/d.tsv", "C": "/tmp/c.tsv"}, state, events)
    assert events == [("attach", "D", id(state["nodes"]))]     # same grove, same cache, D only
    assert state["cohorts"] == {"C", "D"} and g["COHORTS"] == ["C", "D"]

    events.clear()
    _run_preamble({"C": "/tmp/c.tsv"}, state, events)
    assert events == []                                          # all present: pure memo hit

    events.clear()
    g = _run_preamble({"C": "/tmp/c.tsv"}, state, events, {"BRCA-US": "/tmp/brca.tsv"})
    assert events == [("sv", "BRCA-US")]                         # an SV cohort joins the same grove
    assert state["sv_cohorts"] == {"BRCA-US"} and g["SV_COHORTS"] == ["BRCA-US"]

    events.clear()
    state["gg"] = "/tmp/other.gg"                                # a different grove path
    _run_preamble({"C": "/tmp/c.tsv"}, state, events)
    assert events[:2] == ["clear", "deser"]                      # stale entry dropped first


def test_preamble_hides_the_worker_state_and_hands_out_a_read_only_view():
    """The preamble removes its scratch names from the namespace and binds `GROVE` to a view
    that forwards query methods only — so generated code can neither swap the memoised grove
    nor insert into it and change the next question's answer."""
    state, events = {}, []
    g = _run_preamble({"C": "/tmp/c.tsv"}, state, events)
    assert not {"_CANOPY_STATE", "_state", "_grove", "_c", "_links", "_readonly", "_n"} & g.keys()
    view = g["GROVE"]
    assert view.size() == 1 and len(view) == 1               # reads forward
    with pytest.raises(AttributeError):
        view.insert("chr1", None, {})                        # mutators do not
    assert state["grove"].size() == 1                        # ...so the memoised grove is untouched
