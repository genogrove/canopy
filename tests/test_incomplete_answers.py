"""Real execution results must be validated before either adapter interprets them."""
from types import SimpleNamespace

import pytest

from genogrove_canopy import cli, sandbox, serve


@pytest.mark.parametrize("mode", ["run", "worker"])
@pytest.mark.parametrize("stream", ["stdout", "stderr"])
def test_truncated_output_is_an_error_the_model_sees(monkeypatch, tmp_path, mode, stream):
    """A run whose output hit the sandbox cap is reported to the model as an error, never
    handed back as rows: the session keeps no stdout from it, so an answer built on it
    renders nothing — in the CLI and in serve alike."""
    code = "import sys\nprint('10 records:')\n"
    code += "for i in range(10): print(json.dumps({'id': i}), file=sys." + stream + ")\n"
    seen = []

    def fake_answer(question, system_prompt, session, **kw):
        seen.append(session.run(code))
        return "done", []
    monkeypatch.setattr(cli.llm, "answer", fake_answer)
    gg = tmp_path / "unused.gg"
    gg.touch()
    monkeypatch.setattr(serve.resources, "_all_grove_gg", lambda _: gg)
    from contextlib import nullcontext

    with (sandbox.Worker(output_cap=40) if mode == "worker" else nullcontext()) as worker:
        execute = worker.submit if mode == "worker" else lambda s: sandbox.run(s, output_cap=40)
        args = SimpleNamespace(model="test", show_code=False, cohort=None, format="text", trace=None)
        out, error, *_ = cli._answer("question", system_prompt="", base="", gg=str(gg),
                                     args=args, execute=execute)
        assert seen[-1].startswith("error:") and "incomplete" in seen[-1]
        assert out == "done\n" and not error
        grove = SimpleNamespace(system_prompt="", model="test", gg=str(gg), preamble="",
                                worker=SimpleNamespace(submit=execute))
        monkeypatch.setattr(serve, "_grove", lambda _: grove)
        result = serve._pipeline("question", "", "test", lambda *a: None)
        assert "incomplete" in seen[-1]
        assert result["answer"] == "done" and result["records"] == [] and result["summary"] == []


def test_success_and_failure_validation():
    assert sandbox.result_error(sandbox.SandboxResult("answer", "", 0)) == ""
    assert sandbox.result_error(sandbox.SandboxResult("partial", "failure", 7)) == "failure"
    assert sandbox.result_error(sandbox.SandboxResult("partial", "", 7)) == "(the generated code failed with no output)"
    assert sandbox.result_error(sandbox.SandboxResult("partial", "[worker] killed", -1, timed_out=True)) == "[worker] killed"
    # a failure that also overflowed: the traceback stays first, the overflow is a footnote
    both = sandbox.result_error(sandbox.SandboxResult("partial", "Traceback: boom", 1, truncated=True))
    assert both.startswith("Traceback: boom") and "cut at the sandbox cap" in both


def test_failure_after_overflow_shows_the_traceback(tmp_path):
    """Real run: fill stdout past the cap, then raise — the user must see the exception. The cap
    is per stream, so a 1 KB cap cuts the 2.2 KB stdout flood but leaves the traceback whole."""
    r = sandbox.run("for i in range(200): print('x' * 10)\nraise RuntimeError('boom')\n", output_cap=1024)
    assert r.returncode != 0 and r.truncated
    msg = sandbox.result_error(r)
    assert "RuntimeError: boom" in msg and msg.rstrip().endswith("(its output was also cut at the sandbox cap)")
