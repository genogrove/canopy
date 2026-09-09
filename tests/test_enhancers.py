# SPDX-License-Identifier: GPL-3.0-or-later
"""Enhancer layer: declaration parsing, the sandbox preamble (attach-into-grove, memoisation,
the read-only view), `attach_links` against real bindings, and the pinned index bookkeeping."""
import pytest

from genogrove_canopy import llm
from genogrove_canopy import preamble
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


def test_links_file_builds_the_plain_table_the_sandbox_reads(tmp_path, monkeypatch):
    """`links_file` turns a cohort's gzipped byEnhancer index into the plain TSV the sandbox is
    granted: exactly the ten columns `attach_links` reads, the Ensembl id unversioned, the target
    TSS filled in from the bundled gene table, rows whose target is unknown dropped here (the
    sandbox could not resolve them anyway), and the result cached."""
    import gzip

    cohort = "EFO:0000001"
    monkeypatch.setattr(enhancers, "INDEX_DIR", tmp_path / "idx")
    monkeypatch.setattr(enhancers, "LINKS_DIR", tmp_path / "links")
    monkeypatch.setattr(enhancers, "ensure_index", lambda c: True)
    src = enhancers.INDEX_DIR / f"{enhancers._slug(cohort)}.byEnhancer.tsv.gz"
    src.parent.mkdir()
    row = lambda ens, start: "\t".join(  # noqa: E731 — the 11 index columns, in _FIELDS order
        ("chr8", str(start), str(start + 500), "MYC", ens, "intergenic", "FALSE",
         cohort, "2", "0.4", "0.9")) + "\n"
    with gzip.open(src, "wt") as fh:
        fh.write(row("ENSG00000136997.20", 100) + row("ENSG99999999", 200) + row("ENSG00000136997", 300))

    dest = enhancers.links_file(cohort)
    rows = [ln.split("\t") for ln in dest.read_text().splitlines()]
    assert dest == enhancers.LINKS_DIR / f"{enhancers._slug(cohort)}.links.tsv"
    assert [r[1] for r in rows] == ["100", "300"]                    # the unknown target is dropped
    assert rows[0] == ["chr8", "100", "600", "intergenic", "ENSG00000136997",
                       "chr8", "127735434", "2", "0.4", "0.9"]         # unversioned id, TSS filled in
    assert not list(enhancers.LINKS_DIR.glob("*.tmp"))               # written atomically

    src.unlink()
    assert enhancers.links_file(cohort) == dest                      # cached: the index is not reread
