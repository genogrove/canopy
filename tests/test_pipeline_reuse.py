"""Exercise production preamble selection through both warm-worker adapters."""
from types import SimpleNamespace

import pytest

pg = pytest.importorskip("pygenogrove")
from genogrove_canopy import cli, preamble, sandbox, serve
from genogrove_canopy.layers import sv


@pytest.mark.parametrize("adapter", ["cli", "web"])
def test_layer_plain_layer_questions_reuse_one_attached_grove(tmp_path, monkeypatch, adapter):
    gg = (tmp_path / "mini.gg").resolve()
    g = pg.Grove()
    g.insert("chr1", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "A"})
    g.serialize(str(gg))
    tables = {}
    for cohort, sample in [("BRCA-US", "S1"), ("BRCA-UK", "S2")]:
        table = (tmp_path / (cohort + ".tsv")).resolve()
        table.write_text("\t".join(sv._FIELDS + ("sample", "cohort")) + "\n" +
                         f"chr1\t1500\t1501\tchr1\t9000000\t9000001\tSV1\t5\t+\t-\tDEL\tm\t{sample}\t{cohort}\n")
        tables[cohort] = table
    monkeypatch.setattr(sv, "cohort_file", tables.__getitem__)
    monkeypatch.setattr(serve.resources, "_all_grove_gg", lambda _: gg)
    probe = ("gene = next(iter(GROVE.intersect(pg.GenomicCoordinate('*', 1500, 1500), 'chr1')))\n"
             "print(json.dumps({'samples': sorted(m['sample'] for m in GROVE.get_edges(gene) "
             "if m and m.get('rel') == 'breakpoint_edge')}))\n")
    def fake_answer(question, system_prompt, session, **kw):
        if question:
            session.attach(["sv"], question)
            session.run(probe)
        else:
            session.run("print('plain')\n")
        return "answer", []
    monkeypatch.setattr(cli.llm, "answer", fake_answer)
    base = preamble.build(str(gg))
    with sandbox.Worker(data_paths=[tmp_path], extra_syspath=[cli._pygenogrove_site_dir()]) as worker:
        monkeypatch.setattr(serve, "_grove", lambda _: SimpleNamespace(
            gg=str(gg), preamble=base, system_prompt="", model="test", worker=worker))
        args = SimpleNamespace(model="test", show_code=False, cohort=None, format="json", trace=None)
        def ask(question):
            if adapter == "web":
                result = serve._pipeline(question, "", "test", lambda *a: None)
                assert "error" not in result, result
                return result.get("records")
            out, err, *_ = cli._answer(question, system_prompt="", base=base, gg=str(gg),
                                       args=args, execute=worker.submit)
            assert not err, err
            return cli._parse_output(out)[0]
        assert ask("BRCA-US") == [{"samples": ["S1"]}]
        assert ask("") == []                                  # a plain question in between must not evict the memo
        assert ask("BRCA-UK") == [{"samples": ["S1", "S2"]}]
        assert ask("BRCA-US") == [{"samples": ["S1", "S2"]}]


def test_build_puts_own_data_into_the_warm_grove_for_later_runs(tmp_path):
    """The agent's `build` reads a granted own file, inserts into the memoised grove, and a plain
    `run` of the same session — and of the next question's session — sees the record through the
    read-only view; a `run` cannot insert. Real worker, real bindings."""
    gg = (tmp_path / "mini.gg").resolve()
    g = pg.Grove()
    g.insert("chr1", pg.GenomicCoordinate("+", 1000, 2000), {"type": "gene", "id": "A"})
    g.serialize(str(gg))
    mine = (tmp_path / "mine.bed").resolve()
    mine.write_text("chr1\t1400\t1600\tpeak1\n")
    build = (f"for ln in open({str(mine)!r}):\n"
             "    c, s, e, n = ln.split()\n"
             "    k = GROVE.insert(c, pg.GenomicCoordinate('*', int(s), int(e) - 1), {'type': 'peak', 'name': n, 'source': 'own'})\n"
             "    for gk in GROVE.intersect(pg.GenomicCoordinate('*', int(s), int(e) - 1), c):\n"
             "        if gk.data.get('type') == 'gene': GROVE.add_edge(k, gk, {'rel': 'overlaps', 'source': 'own'})\n"
             "print('built: 1')\n")
    probe = ("hits = [k.data for k in GROVE.intersect(pg.GenomicCoordinate('*', 1500, 1500), 'chr1')]\n"
             "print(json.dumps(sorted(d.get('name') or d['id'] for d in hits)))\n")
    base = preamble.build(str(gg))
    with sandbox.Worker(data_paths=[gg, mine], extra_syspath=[cli._pygenogrove_site_dir()]) as worker:
        s1 = cli.Session(gg=str(gg), base=base, execute=worker.submit, say=lambda m: None)
        assert s1.run(probe) == '["A"]\n'                      # before the build: lazy view, gene only
        assert s1.build(build) == "built: 1\n"
        assert s1.run(probe) == '["A", "peak1"]\n'             # same session: warm view sees the peak
        assert s1.run("GROVE.insert('chr1', pg.GenomicCoordinate('*', 1, 2), {})\n").startswith("error:")
        # a build may add, never take away: removal and clearing are off the handle
        err = s1.build("GROVE.clear_graph()\n")
        assert err.startswith("error:") and "build allows insert, add_edge, add_external_key" in err
        assert s1.build("GROVE.remove_key(next(iter(GROVE.intersect(pg.GenomicCoordinate('*', 1500, 1500), 'chr1'))))\n").startswith("error:")
        assert s1.run(probe) == '["A", "peak1"]\n'
        # the next question's session: with own files given, it starts warm — the memo kept the peak
        s2 = cli.Session(gg=str(gg), base=base, execute=worker.submit, say=lambda m: None, warm=True)
        assert s2.run(probe) == '["A", "peak1"]\n'
