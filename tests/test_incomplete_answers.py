"""Real execution results must be validated before either adapter renders them."""
from types import SimpleNamespace

import pytest

from genogrove_canopy import cli, sandbox, serve


@pytest.mark.parametrize("mode", ["run", "worker"])
@pytest.mark.parametrize("stream", ["stdout", "stderr"])
def test_truncated_output_is_an_error_in_both_adapters(monkeypatch, tmp_path, mode, stream):
    code = "import sys\nprint('10 records:')\n"
    code += "for i in range(10): print(json.dumps({'id': i}), file=sys." + stream + ")\n"
    monkeypatch.setattr(cli.llm, "generate_query", lambda *a, **kw: ("", [], code))
    gg = tmp_path / "unused.gg"
    gg.touch()
    monkeypatch.setattr(serve.resources, "_all_grove_gg", lambda _: gg)
    with sandbox.Worker(output_cap=40) as worker:
        execute = worker.submit if mode == "worker" else lambda s: sandbox.run(s, output_cap=40)
        args = SimpleNamespace(model="test", show_code=False, cohort=None, format="json")
        out, error, *_ = cli._answer("question", system_prompt="", base="", gg=str(gg),
                                     args=args, execute=execute)
        assert not out and "incomplete" in error
        grove = SimpleNamespace(system_prompt="", model="test", gg=str(gg), preamble="",
                                worker=SimpleNamespace(submit=execute))
        monkeypatch.setattr(serve, "_grove", lambda _: grove)
        result = serve._pipeline("question", "", "test", lambda *a: None)
        assert "incomplete" in result["error"]
        assert "summary" not in result and "records" not in result


def test_success_and_failure_validation():
    assert sandbox.result_error(sandbox.SandboxResult("answer", "", 0)) == ""
    assert sandbox.result_error(sandbox.SandboxResult("partial", "failure", 7)) == "failure"
    assert "time limit" in sandbox.result_error(sandbox.SandboxResult("partial", "", -1, timed_out=True))
