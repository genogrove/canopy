# SPDX-License-Identifier: GPL-3.0-or-later
"""The sandbox must actually run pygenogrove: import it (despite -S) and
deserialize a .gg. Validates the extra_syspath fix + the host/agent .gg path.
Runs only where pygenogrove is installed (CI skips)."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from genogrove_canopy import sandbox
from genogrove_canopy.cli import _pygenogrove_site_dir
from genogrove_canopy.gff import load_gff

GFF3 = (
    "##gff-version 3\n"
    "chr1\tH\tgene\t1000\t2000\t.\t+\t.\tID=g1;gene_name=AAA\n"
    "chr1\tH\ttranscript\t1000\t2000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tH\texon\t1000\t1500\t.\t+\t.\tID=e1;Parent=t1\n"
).encode()


def test_sandbox_deserializes_and_queries_a_grove(tmp_path):
    src = tmp_path / "mini.gff3"
    src.write_bytes(GFF3)
    gg = tmp_path / "mini.gg"
    load_gff(src).serialize(str(gg))

    # This is the shape of what the CLI runs: a host-injected path var + agent code.
    code = (
        f"GG = {str(gg)!r}\n"
        "import pygenogrove as pg\n"
        "g = pg.Grove.deserialize(GG)\n"
        "q = pg.GenomicCoordinate('*', 1200, 1200)\n"
        "print('genes', sum(1 for k in g.intersect(q, 'chr1') if k.data['type'] == 'gene'))\n"
        "print('size', g.size())\n"
    )
    result = sandbox.run(
        code, data_paths={"g": str(gg)}, extra_syspath=[_pygenogrove_site_dir()]
    )

    assert result.returncode == 0, result.stderr
    assert "genes 1" in result.stdout
    assert "size 1" in result.stdout  # only the gene is indexed; transcript + exon external


def test_sandbox_still_blocks_network_with_pygenogrove_on_path(tmp_path):
    # The widened sys.path must not let the untrusted code import a blocked module.
    result = sandbox.run("import socket", extra_syspath=[_pygenogrove_site_dir()])
    assert result.returncode != 0
    assert "blocked" in result.stderr or "socket" in result.stderr


def test_enhancer_preamble_attaches_links_inside_the_sandbox(tmp_path, monkeypatch):
    """The cohort preamble opens a links table from `enhancers.LINKS_DIR` inside the sandbox,
    so `_grove_context` must grant that directory — with only the grove granted, every enhancer
    question died on `PermissionError` and no test executed the preamble to notice."""
    from genogrove_canopy import resources
    from genogrove_canopy.cli import _grove_context
    from genogrove_canopy.layers import enhancers

    gg = tmp_path / "mini.gg"
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 999, 1999), {"type": "gene", "id": "ENSG1.2"})
    g.serialize(str(gg))
    links = tmp_path / "links" / "c.links.tsv"
    links.parent.mkdir()
    links.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t1\t0.5\t0.9\n")
    monkeypatch.setattr(enhancers, "LINKS_DIR", links.parent)
    monkeypatch.setattr(resources, "ensure_all_grove", lambda name: gg)

    _, _, data_paths = _grove_context()
    code = enhancers.preamble(str(gg), {"COH": str(links)}) + (
        "enh = [k for k in GROVE.intersect(pg.GenomicCoordinate('*', 150, 150), 'chr1')]\n"
        "print('enh', len(enh), enh[0].data['source'])\n"
    )
    result = sandbox.run(code, data_paths=data_paths, extra_syspath=[_pygenogrove_site_dir()])

    assert result.returncode == 0, result.stderr
    assert "enh 1 ENCODE-rE2G" in result.stdout


def test_warm_worker_reuses_the_attached_grove_and_a_query_cannot_mutate_it(tmp_path):
    """The path `serve` and `-i` always take: the cohort preamble runs twice in one `Worker`.
    The second run must hit `_CANOPY_STATE` (no re-attach) and see exactly what the first saw —
    even though the first query tried to insert into `GROVE`, which the read-only view refuses."""
    from genogrove_canopy.layers import enhancers

    gg = tmp_path / "mini.gg"
    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 999, 1999), {"type": "gene", "id": "ENSG1.2"})
    g.serialize(str(gg))
    links = tmp_path / "c.links.tsv"
    links.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t1\t0.5\t0.9\n")
    pre = "import json\n" + enhancers.preamble(str(gg), {"C": str(links)})
    count = ("enh = [k for k in GROVE.intersect(pg.GenomicCoordinate('*', 0, 5000), 'chr1')"
             " if k.data.get('type') == 'enhancer']\nprint('enh', len(enh))\n")

    # The view hides the grove, but a bound method's `__self__` still names it — good enough to
    # prove both runs used the same object, i.e. the second was a memo hit.
    ident = "print('id', id(GROVE.intersect.__self__))\n"
    links2 = tmp_path / "d.links.tsv"       # a second cohort: the same element plus a new one
    links2.write_text("chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t2\t0.4\t0.7\n"
                      "chr1\t300\t400\tintergenic\tENSG1\tchr1\t1000\t1\t0.1\t0.2\n")
    pre2 = "import json\n" + enhancers.preamble(str(gg), {"D": str(links2)})
    both = ("gene = next(iter(GROVE.intersect(pg.GenomicCoordinate('*', 1500, 1500), 'chr1')))\n"
            "print('cohorts', sorted(set(c for m in GROVE.get_edges(gene)"
            " if m and m.get('rel') == 'regulated_by' for c in m['byCohort'])), COHORTS)\n")
    w = sandbox.Worker(data_paths=[str(tmp_path)], extra_syspath=[_pygenogrove_site_dir()])
    try:
        first = w.submit(pre + count + ident
                         + "try:\n    GROVE.insert('chr1', pg.GenomicCoordinate('.', 5, 6), {})\n"
                         "except AttributeError as e:\n    print('refused', e)\n")
        # A plain (no-cohort) question in between must not be able to reach the memo either.
        plain = w.submit("import json\n" + enhancers.preamble(str(gg))
                         + "print('hidden', '_CANOPY_STATE' not in globals())\n")
        second = w.submit(pre + count + ident)
        third = w.submit(pre2 + count + ident + both)
    finally:
        w.close()

    assert first.returncode == 0, first.stderr
    assert "enh 1" in first.stdout and "refused" in first.stdout
    assert plain.returncode == 0 and "hidden True" in plain.stdout, plain.stderr
    assert second.returncode == 0, second.stderr
    assert "enh 1" in second.stdout                       # the refused insert left no trace
    assert third.returncode == 0, third.stderr
    # Cohort D attached onto the SAME grove: its new element appears, the shared element merged
    # (one node, both cohorts on the gene's edges), and COHORTS names only this question's.
    assert "enh 2" in third.stdout
    assert "cohorts ['C', 'D'] ['D']" in third.stdout
    ids = [ln for ln in (first.stdout + second.stdout + third.stdout).splitlines() if ln.startswith("id ")]
    assert len(ids) == 3 and len(set(ids)) == 1            # same grove object throughout


def test_system_prompt_worked_example_runs_against_the_attached_grove(tmp_path):
    """The prompt's worked example is the codegen contract in executable form. Run it, verbatim,
    through the real sandbox on a small grove with a cohort attached — the drift this guards
    against (prompt says one thing, preamble does another) shipped enhancer answers with zero
    rows once, and no test executed the example."""
    import json
    import re
    from pathlib import Path

    from genogrove_canopy.layers import enhancers

    md = (Path(enhancers.__file__).parents[1] / "prompts" / "system.md").read_text()
    m = re.search(r"### Worked example.*?COHORT: MCF-7\nTARGETS: .*?\n\n```python\n(.*?)```", md, re.S)
    assert m, "worked example not found in system.md"
    example = m.group(1)

    # EGFR at the example's variant (chr7:55,191,822 -> closed 55_191_821), one cCRE under the
    # enhancer window, one MCF-7 link whose target TSS (1-based) lies inside the gene.
    g = pg.Grove(order=100)
    g.insert("chr7", pg.GenomicCoordinate("+", 55_018_819, 55_211_627),
             {"type": "gene", "id": "ENSG00000146648.23", "name": "EGFR", "biotype": "protein_coding"})
    g.insert("chr7", pg.GenomicCoordinate(".", 55_018_400, 55_018_700),
             {"type": "regulatory_region", "source": "ENCODE-SCREEN", "class": "PLS",
              "id": "EH38E0000001", "rdhs": "EH38D0000001"})
    gg = tmp_path / "egfr.gg"
    g.serialize(str(gg))
    links = tmp_path / "mcf7.links.tsv"
    links.write_text("chr7\t55018383\t55019773\tpromoter\tENSG00000146648\tchr7\t55018820\t1\t0.9\t0.99999\n")

    code = "import json\n" + enhancers.preamble(str(gg), {"EFO:0001203": str(links)}) + example
    result = sandbox.run(code, data_paths=[str(tmp_path)], extra_syspath=[_pygenogrove_site_dir()])

    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert lines[0] == "variant chr7:55,191,822 in EGFR (1 gene, 1 enhancer links):"
    rows = [json.loads(ln) for ln in lines[1:]]
    assert [r["type"] for r in rows] == ["gene", "enhancer"]
    enh = rows[1]
    assert (enh["start"], enh["end"]) == (55_018_383, 55_019_772)       # BED half-open -> closed
    assert enh["score"] == 0.99999 and enh["n"] == 1 and enh["cohort"] == "EFO:0001203"
    assert enh["target"] == "EGFR" and enh["name"] == "enh:promoter->EGFR"
    assert enh["ccre_overlap"] == [{"id": "EH38E0000001", "class": "PLS", "bp": 301}]


def test_grove_context_grants_resolved_paths_so_a_symlinked_cache_works(tmp_path, monkeypatch):
    """The sandbox resolves its granted roots; if the host handed the preamble a path spelled
    through a symlink (a `GENOGROVE_CANOPY_CACHE` under `/var`, say), every read was refused."""
    from genogrove_canopy import resources
    from genogrove_canopy.cli import _grove_context
    from genogrove_canopy.layers import enhancers

    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real, target_is_directory=True)

    g = pg.Grove(order=100)
    g.insert("chr1", pg.GenomicCoordinate("+", 999, 1999), {"type": "gene", "id": "ENSG1.2"})
    g.serialize(str(real / "mini.gg"))
    (real / "links").mkdir()
    (real / "links" / f"{enhancers._slug('C')}.links.tsv").write_text(  # the name links_file expects
        "chr1\t100\t200\tintergenic\tENSG1\tchr1\t1000\t1\t0.5\t0.9\n")
    monkeypatch.setattr(resources, "ensure_all_grove", lambda name: link / "mini.gg")   # via the symlink
    monkeypatch.setattr(enhancers, "LINKS_DIR", link / "links")
    monkeypatch.setattr(enhancers, "ensure_index", lambda c: True)

    _, base_pre, data_paths = _grove_context()
    links = enhancers.links_file("C")
    assert "link/" not in base_pre and all("link/" not in p for p in data_paths) and "link/" not in str(links)

    code = (base_pre + enhancers.preamble(data_paths[0], {"C": str(links)})
            + "print('enh', sum(1 for k in GROVE.intersect(pg.GenomicCoordinate('*', 150, 150), 'chr1')"
              " if k.data.get('type') == 'enhancer'))\n")
    result = sandbox.run(code, data_paths=data_paths, extra_syspath=[_pygenogrove_site_dir()])
    assert result.returncode == 0, result.stderr
    assert "enh 1" in result.stdout
