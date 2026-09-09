"""Filesystem regressions through ordinary Python and native APIs, in both modes."""
from pathlib import Path

import pytest

from genogrove_canopy import sandbox


@pytest.fixture(params=["run", "worker"])
def execute(request, tmp_path):
    allowed = (tmp_path / "allowed").resolve()
    allowed.mkdir()
    options = {"data_paths": [allowed]}
    if request.param == "run":
        yield lambda code: sandbox.run(code, **options), allowed
    else:
        with sandbox.Worker(**options) as worker:
            yield worker.submit, allowed


@pytest.mark.parametrize("api", ["open", "io.open", "io.FileIO"])
def test_python_io_cannot_read_or_truncate_ungranted_file(execute, api):
    run, allowed = execute
    sentinel = allowed.parent / "sentinel"
    sentinel.write_text("KEEP")
    for mode in ["r", "w"]:
        result = run(f"import io\n{api}({str(sentinel)!r}, {mode!r}).close()")
        assert result.returncode != 0
        assert sentinel.read_text() == "KEEP"
    data = allowed / "data"
    data.write_text("READ")
    result = run(f"import io\nprint({api}({str(data)!r}, 'r').read())")
    assert result.returncode == 0 and "READ" in result.stdout, result
    result = run(f"import io\n{api}({str(data)!r}, 'w').close()")
    assert result.returncode != 0 and data.read_text() == "READ"


def test_symlink_does_not_extend_granted_directory(execute):
    run, allowed = execute
    sentinel = allowed.parent / "sentinel"
    sentinel.write_text("KEEP")
    link = allowed / "link"
    link.symlink_to(sentinel)
    result = run(f"import io\nprint(io.open({str(link)!r}).read())")
    assert result.returncode != 0 and "KEEP" not in result.stdout


@pytest.mark.parametrize("mode", ["run", "worker"])
def test_native_reader_respects_grants_and_symlinks(tmp_path, mode):
    pg = pytest.importorskip("pygenogrove")
    from genogrove_canopy.cli import _pygenogrove_site_dir
    allowed = (tmp_path / "allowed").resolve()
    allowed.mkdir()
    outside = (tmp_path / "outside.gff3").resolve()
    outside.write_text("##gff-version 3\nchr1\tH\tgene\t1\t10\t.\t+\t.\tID=sentinel\n")
    inside = allowed / "inside.gff3"
    inside.write_text(outside.read_text())
    link = allowed / "link.gff3"
    link.symlink_to(outside)
    options = dict(data_paths=[allowed], extra_syspath=[_pygenogrove_site_dir()])
    with sandbox.Worker(**options) as worker:
        run = worker.submit if mode == "worker" else lambda s: sandbox.run(s, **options)
        for path, permitted in [(outside, False), (link, False), (inside, True)]:
            result = run(f"import pygenogrove as pg\nprint(len(list(pg.GffReader({str(path)!r}))))")
            assert (result.returncode == 0) == permitted, result
            if permitted:
                assert result.stdout.strip() == "1"
        before = inside.read_bytes()
        result = run(f"import pygenogrove as pg\npg.Grove().serialize({str(inside)!r})")
        assert result.returncode != 0 and inside.read_bytes() == before


def test_unsupported_platform_fails_before_query_runs(tmp_path):
    import subprocess
    import sys
    # Execute the trusted bootstrap with its platform detection replaced in a child only.
    script = sandbox._build_script("print('QUERY RAN')", [], [])
    script = script.replace('if sys.platform == "darwin":', 'if False:').replace(
        'elif sys.platform == "linux" and os.uname().machine in ("x86_64", "aarch64"):', 'elif False:')
    child = subprocess.run([sys.executable, "-I", "-S", "-c", script], capture_output=True, text=True)
    assert child.returncode != 0 and "QUERY RAN" not in child.stdout
    assert "Filesystem isolation requires" in child.stderr


def test_worker_reports_isolation_initialization_failure(monkeypatch):
    monkeypatch.setattr(sandbox.inspect, "getsource", lambda _: (
        "def restrict_filesystem(roots):\n"
        "    raise RuntimeError('test isolation unavailable')\n"))
    with pytest.raises(RuntimeError, match="test isolation unavailable"):
        sandbox.Worker()
    result = sandbox.run("print('QUERY RAN')")
    assert result.returncode != 0 and "QUERY RAN" not in result.stdout
    assert "test isolation unavailable" in result.stderr
