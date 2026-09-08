# SPDX-License-Identifier: GPL-3.0-or-later
"""Host-side enhancer layer: target parsing, gene→TSS resolution, the ENHANCERS injection,
and (when the index bundle is present) the tabix lookups. The generated query code can't read
the index — the sandbox allowlist blocks file I/O — so these run host-side only."""
import pytest

from genogrove_canopy import llm
from genogrove_canopy.layers import enhancers

FLAGSHIP = "EFO:0005726"  # LNCaP


def test_parse_cohort_targets_and_code():
    text = ('reasoning...\nCOHORT: MCF-7\n'
            'TARGETS: [{"gene": "MYC"}, {"region": "chr8:127700000-127740000"}]\n'
            "```python\nprint(1)\n```")
    cohort, targets, code = llm.parse_targets_and_code(text)
    assert cohort == "MCF-7"
    assert targets == [{"gene": "MYC"}, {"region": "chr8:127700000-127740000"}]
    assert code == "print(1)\n"


def test_parse_none_declared_is_empty():
    # a structural (non-enhancer) reply declares nothing -> backward compatible
    cohort, targets, code = llm.parse_targets_and_code("```python\nx = 1\n```")
    assert cohort == "" and targets == [] and code == "x = 1\n"


def test_parse_malformed_targets_tolerated():
    cohort, targets, code = llm.parse_targets_and_code("TARGETS: [not json\n```python\np()\n```")
    assert cohort == "" and targets == [] and code == "p()\n"


@pytest.mark.parametrize("text", [
    # declarations in a BARE fence, program in a python fence (the real crash: 'MCF-7' -> NameError)
    '```\nCOHORT: MCF-7\nTARGETS: [{"gene": "EGFR"}]\n```\n```python\np()\n```',
    # declarations leaked INSIDE the python fence as leading lines
    '```python\nCOHORT: MCF-7\nTARGETS: [{"gene": "EGFR"}]\np()\n```',
    # clean: plain declaration lines + one python fence
    'COHORT: MCF-7\nTARGETS: [{"gene": "EGFR"}]\n```python\np()\n```',
])
def test_declarations_never_leak_into_code(text):
    cohort, targets, code = llm.parse_targets_and_code(text)
    assert cohort == "MCF-7" and targets == [{"gene": "EGFR"}]
    assert "COHORT" not in code and "TARGETS" not in code and code.strip() == "p()"


def test_resolve_gene():
    # gene_tss.tsv.gz ships in the package, so this works without the index
    hit = enhancers.resolve_gene("MYC")
    assert hit == ("ENSG00000136997", "chr8", 127735434)
    assert enhancers.resolve_gene("NOT_A_GENE") is None


def test_preamble_no_cohorts_opens_lazily():
    pre = enhancers.preamble("/tmp/x.gg")
    assert pre == ('import pygenogrove as pg\n'
                    'GROVE = pg.GroveView.open("/tmp/x.gg")\nCOHORTS = []\nENHANCERS = []\n')
    compile(pre, "<preamble>", "exec")


def test_preamble_with_cohorts_attaches_each_onto_one_grove():
    pre = enhancers.preamble(
        "/tmp/x.gg", {"EFO:0005726": "/tmp/a.tsv", "EFO:0009318": "/tmp/b.tsv"})
    assert 'grove=pg.Grove.deserialize("/tmp/x.gg")' in pre
    assert '[["EFO:0005726", "/tmp/a.tsv"], ["EFO:0009318", "/tmp/b.tsv"]]' in pre
    assert "attach_links(_grove, _links, _c, _state['nodes'])" in pre
    assert 'COHORTS = ["EFO:0005726", "EFO:0009318"]' in pre
    assert "GROVE = _readonly(_grove)" in pre
    assert "_CANOPY_STATE" in pre
    compile(pre, "<preamble>", "exec")


@pytest.mark.skipif(not enhancers.index_present(FLAGSHIP),
                    reason="rE2G index bundle not present (download/build it first)")
def test_fetch_for_targets_gene_and_region():
    recs = enhancers.fetch_for_targets(
        [{"gene": "MYC"}, {"region": "chr8:127700000-127740000"}], [FLAGSHIP])
    assert recs
    assert all("target_gene" in r and "score_max" in r for r in recs)
    # deduped across the two overlapping targets
    keys = {(r["chrom"], r["start"], r["end"], r["target_gene"], r["cohort"]) for r in recs}
    assert len(keys) == len(recs)


def test_index_present_requires_every_file_not_just_the_tables(tmp_path, monkeypatch):
    """A cohort with tables but no tabix indexes is *not* ready.

    `index_present` used to check only the two `.tsv.gz` files, so an interrupted first fetch —
    or a cache from the old locally-built layout — reported itself ready and then failed inside
    tabix. `_index_files` is now the single definition of the set, used by the check and the
    fetch alike.
    """
    monkeypatch.setattr(enhancers, "INDEX_DIR", tmp_path)
    names = enhancers._index_files("EFO:0005726")
    assert len(names) == 4, "a cohort needs two tables and their two indexes"

    for name in names:
        if name.endswith(".tbi"):
            continue
        (tmp_path / name).write_bytes(b"")
    assert not enhancers.index_present("EFO:0005726"), "tables alone must not count as ready"

    for name in names:
        (tmp_path / name).write_bytes(b"")
    assert enhancers.index_present("EFO:0005726")


def test_attach_links_merges_a_second_cohort_onto_one_node_and_edge(tmp_path):
    """The same element→gene link in two cohorts is one enhancer node and one `regulates` edge
    whose `byCohort` carries both — not a second node at the same interval with a parallel edge,
    which is what two plain `add_edge` calls (and a per-call node cache) produced."""
    pg = pytest.importorskip("pygenogrove")

    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 999, 1999), {"type": "gene", "id": "ENSG1.2"})
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    a.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t1\t0.5\t0.9\n")
    b.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t2\t0.4\t0.7\n"
                 "chr1\t300\t400\tintergenic\tENSG1\tchr1\t1000\t1\t0.1\t0.2\n")
    nodes = {}  # shared across cohorts, as the preamble does
    assert enhancers.attach_links(g, a, "C1", nodes) == (1, 1, 0)
    assert enhancers.attach_links(g, b, "C2", nodes) == (2, 2, 0)

    gene = next(k for k in g.intersect(pg.GenomicCoordinate("*", 1500, 1500), "chr1"))
    enh = [k for k in g.intersect(pg.GenomicCoordinate("*", 150, 150), "chr1")
           if k.data.get("type") == "enhancer"]
    assert len(enh) == 1
    edges = [m for _, m in g.get_edge_list(enh[0])]
    assert edges == [{"rel": "regulates", "byCohort": {
        "C1": {"score_max": 0.9, "score_mean": 0.5, "n_rep": 1},
        "C2": {"score_max": 0.7, "score_mean": 0.4, "n_rep": 2}}}]
    back = [m for t, m in g.get_edge_list(gene) if t.value.start == 100]
    assert len(back) == 1 and back[0]["byCohort"].keys() == {"C1", "C2"}
    assert len(g.get_neighbors_if(gene, lambda m: m and m.get("rel") == "regulated_by")) == 2


def _run_preamble(cohort_links, state, events):
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
    pre = enhancers.preamble("/tmp/x.gg", cohort_links)
    body = pre.split("import pygenogrove as pg\n", 1)[1].replace(
        "        attach_links(_grove, _links, _c, _state['nodes'])",
        "        events.append(('attach', _c, id(_state['nodes'])))")
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
