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
    monkeypatch.setattr(cli.llm, "generate_query", lambda q, *a, **kw:
                        (q, ["sv"], probe) if q else ("", [], "print('plain')\n"))
    base = preamble.build(str(gg))
    with sandbox.Worker(data_paths=[tmp_path], extra_syspath=[cli._pygenogrove_site_dir()]) as worker:
        monkeypatch.setattr(serve, "_grove", lambda _: SimpleNamespace(
            gg=str(gg), preamble=base, system_prompt="", model="test", worker=worker))
        args = SimpleNamespace(model="test", show_code=False, cohort=None, format="json")
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
        ask("")
        assert ask("BRCA-UK") == [{"samples": ["S1", "S2"]}]
        assert ask("BRCA-US") == [{"samples": ["S1", "S2"]}]
