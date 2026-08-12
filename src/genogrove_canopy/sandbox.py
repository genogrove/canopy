# SPDX-License-Identifier: GPL-3.0-or-later
"""Sandboxed execution of LLM-generated Python.

**Security-critical.** The code passed to :func:`run` is produced by a language
model and is treated as untrusted. It is executed out-of-process under several
layers of restriction.

Threat model
------------
The realistic adversary here is *LLM-generated code* — code the model wrote to
answer a genomics question, possibly buggy or accidentally dangerous (a runaway
loop, an unintended network call, a stray file write), not a human deliberately
crafting an escape. The guarantees below are sized to that model.

What is enforced by the parent / the OS (a child cannot lift these)
-------------------------------------------------------------------
* **Out-of-process.** The code runs in a separate interpreter (``subprocess``),
  never via in-process ``exec``/``eval``.
* **Stripped environment.** The child gets a minimal env — no secrets such as
  ``ANTHROPIC_API_KEY`` are reachable even if the in-child guards are defeated.
* **Resource caps** (POSIX ``setrlimit`` in a pre-exec hook): CPU seconds,
  address space (memory), ``RLIMIT_FSIZE = 0`` (no file *writes* of any size),
  and an open-file-descriptor cap. These can only be *lowered* by the child.
* **Wall-clock kill.** The child runs in its own session; on timeout the whole
  process group is killed.
* **Output-size cap.** stdout/stderr are read with a byte cap so a flood cannot
  exhaust the parent's memory.

In-child guards (strong against generated code; defense-in-depth)
-----------------------------------------------------------------
A bootstrap prelude prepended to the code:

* installs an **import allowlist** (``pygenogrove`` + a small compute-only set),
* **scrubs** the dangerous primitives the interpreter preloads (``posix``,
  ``marshal``, and the network/exec extension modules) from ``sys.modules`` so a
  later ``import os`` / ``socket`` / ``subprocess`` is refused at ``find_spec``,
  before any loader runs — so **network has no path through the import system**
  (no ``socket``/``_socket``/``ctypes`` can be imported), and
* replaces ``open`` with a **read-only** variant restricted to registry-resolved
  data roots.

Residual risk
-------------
The in-child guards run in the same interpreter as the untrusted code, so they
are *not* adversary-proof: code that removes the guard from ``sys.meta_path`` and
imports the built-in ``posix``, or reaches the import machinery's private ``os``
reference, can call ``posix.system`` and from there shell out (which is also a
network path). That is out of scope for the in-child layer by design — the
**parent/OS layer** (stripped env, rlimits, ``RLIMIT_FSIZE=0``, session-kill,
output cap) is the hard boundary, and it holds regardless. True isolation against
a hostile child needs an OS-level backend (seccomp-bpf blocking ``execve`` /
``socket`` / network + mount namespaces / a container / an unprivileged jail);
that is the documented next step, and the architecture keeps it pluggable (it
would wrap the same ``subprocess`` invocation). The threat model the in-child
layer is sized to is LLM-generated code, not a human crafting escapes.
"""

from __future__ import annotations

import json
import os
import select
import shutil
import struct
import subprocess
import sys
import tempfile
import threading
import time
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

try:  # POSIX-only; resource caps are skipped (with a weaker guarantee) elsewhere.
    import resource
except ImportError:  # pragma: no cover - non-POSIX (e.g. Windows)
    resource = None  # type: ignore[assignment]


# Compute-only stdlib the generated code may use. Deliberately excludes anything
# that does I/O, networking, subprocess, or dynamic import (os, io helpers,
# socket, subprocess, ctypes, gzip/csv/pickle, importlib, ...). pygenogrove does
# its own file reading, so file/codec modules are intentionally absent. Verified
# (see tests) not to transitively pull a network/exec primitive back in.
_COMPUTE_MODULES = (
    "math",
    "itertools",
    "collections",
    "functools",
    "operator",
    "heapq",
    "bisect",
    "re",
    "json",
)

#: Top-level modules the generated code is permitted to import.
ALLOWED_IMPORTS = frozenset(("pygenogrove",) + _COMPUTE_MODULES)

#: Default caps.
DEFAULT_TIMEOUT_S = 30.0
DEFAULT_OUTPUT_CAP = 256 * 1024  # bytes, per stream
# 4 GiB address space. A mutable `Grove.deserialize` of the shipped grove peaks at ~2.0 GB, so the
# old 2 GiB cap was exactly on the line — fine on macOS (RLIMIT_AS is unsettable there, so the cap
# silently never applied) and a coin-flip OOM on Linux, where it does.
_DEFAULT_MEMORY_BYTES = 4 * 1024 * 1024 * 1024
_MAX_OPEN_FDS = 64


@dataclass
class SandboxResult:
    """Outcome of running generated code in the sandbox."""

    stdout: str
    stderr: str
    returncode: int
    timed_out: bool = False
    #: True if either stream hit :data:`DEFAULT_OUTPUT_CAP` and was truncated.
    truncated: bool = False


# --------------------------------------------------------------------------- #
# The in-child bootstrap. Prepended to the untrusted code; runs first and
# installs the import allowlist / read-only open before the code executes.
# --------------------------------------------------------------------------- #

_BOOTSTRAP = '''\
import sys as _sys, builtins as _builtins

_ALLOW = set(_ALLOW_JSON)
_ROOTS = tuple(_ROOTS_JSON)

# 0) Make allowlisted *installed* packages importable. The child runs with -S
#    (no site processing), so site-packages is not on sys.path; the parent passes
#    the interpreter's site dir(s) here so e.g. `pygenogrove` can load. The import
#    guard installed below still refuses every non-allowlisted top-level import,
#    so widening sys.path does not widen what the untrusted code can import.
_sys.path[:0] = [p for p in _SYSPATH_JSON if p not in _sys.path]

# 1) Pre-warm the allowlist so no later import needs the file-import machinery
#    (which we are about to disable). A module that is not installed is skipped;
#    importing it from user code then fails cleanly via the normal import error.
for _m in _ALLOW:
    try:
        __import__(_m)
    except Exception:
        pass

# 2) Scrub the dangerous primitives (preloaded, or pulled in while warming a
#    trusted module). Popping a name forces any later `import` of it back through
#    the allowlist guard below, which refuses it — so `import os` / `socket` /
#    `subprocess` all fail at the find_spec stage, before any loader runs.
#    (The import machinery keeps its own private os reference for file-based
#    loads; reaching posix.system through it is the documented residual risk —
#    see the module docstring — and is closed only by the OS-level backend.)
for _m in (
    "posix", "marshal", "os", "socket", "_socket",
    "subprocess", "_posixsubprocess", "ctypes", "_ctypes",
):
    _sys.modules.pop(_m, None)

# 3) Import allowlist: anything already cached (interpreter internals + the
#    pre-warmed allowlist) is fine; any *new* top-level import must be allowed.
#    _INFRA is the codec/text machinery the interpreter loads lazily — e.g.
#    text-mode open() pulls in `encodings.<name>` for whatever the locale
#    resolves to (ascii on a bare 3.9 runner). These are pure data transforms
#    (no os/network/subprocess), so allowing them does not widen the boundary.
_INFRA = {"encodings", "codecs", "_codecs"}

class _Guard:
    def find_spec(self, name, path=None, target=None):
        top = name.split(".")[0]
        if name in _sys.modules or top in _ALLOW or top in _INFRA:
            return None
        raise ImportError("import %r is blocked in the sandbox" % name)

_sys.meta_path.insert(0, _Guard())

# 4) Read-only open, restricted to registry-resolved data roots. Writes are
#    refused here (and RLIMIT_FSIZE=0 enforces it at the OS level regardless).
_real_open = _builtins.open

def _norm(p):
    # Canonicalize an absolute path with pure string ops — no os/posixpath (which
    # would re-import os). Collapses "." and ".." so "/root/../etc" cannot escape;
    # a relative path has no leading "/" and so matches no (absolute) root.
    p = str(p)
    parts = []
    for seg in p.split("/"):
        if seg in ("", "."):
            continue
        if seg == "..":
            if parts:
                parts.pop()
            continue
        parts.append(seg)
    return ("/" if p.startswith("/") else "") + "/".join(parts)

def _guarded_open(file, mode="r", *args, **kwargs):
    if any(flag in mode for flag in ("w", "a", "x", "+")):
        raise PermissionError("the sandbox is read-only")
    target = _norm(file)
    if not any(target == r or target.startswith(r + "/") for r in _ROOTS):
        raise PermissionError(
            "reads are restricted to registry data paths: %r" % (file,)
        )
    return _real_open(file, mode, *args, **kwargs)

_builtins.open = _guarded_open

# ----- end bootstrap; untrusted code follows -----
'''


def _build_script(code: str, roots: list[str], syspath: list[str]) -> str:
    """Prepend the bootstrap (with the allowlist, data roots, and site path baked in)."""
    header = (
        f"_ALLOW_JSON = {json.dumps(sorted(ALLOWED_IMPORTS))}\n"
        f"_ROOTS_JSON = {json.dumps(roots)}\n"
        f"_SYSPATH_JSON = {json.dumps(syspath)}\n"
    )
    return header + _BOOTSTRAP + "\n" + code


def _child_env() -> dict[str, str]:
    """A minimal environment for the child — no inherited secrets.

    Only the few variables needed for a Python interpreter to run with sane
    text handling are passed; notably ``ANTHROPIC_API_KEY`` and everything else
    in the parent environment are dropped.
    """
    env = {
        "PATH": "/usr/bin:/bin",
        "LC_ALL": "C.UTF-8",
        "LANG": "C.UTF-8",
        # Hardening flags for the child interpreter.
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
        # Force deterministic UTF-8 text I/O regardless of the host locale (a
        # bare runner can resolve to ascii, which would otherwise need a codec
        # the import guard hasn't pre-loaded). Also aids reproducibility.
        "PYTHONUTF8": "1",
    }
    return env


def _apply_limits(timeout_s: float) -> None:  # pragma: no cover - runs in child
    """Pre-exec hook (POSIX): drop into a new session and cap resources."""
    os.setsid()  # own process group, so the parent can kill the whole tree
    if resource is None:
        return
    cpu = int(timeout_s) + 1
    _setrlimit(resource.RLIMIT_CPU, cpu)
    _setrlimit(resource.RLIMIT_AS, _DEFAULT_MEMORY_BYTES)
    _setrlimit(resource.RLIMIT_FSIZE, 0)  # no file writes of any size
    _setrlimit(resource.RLIMIT_NOFILE, _MAX_OPEN_FDS)


def _setrlimit(which: int, value: int) -> None:  # pragma: no cover - runs in child
    try:
        soft, hard = resource.getrlimit(which)
        new_hard = value if hard == resource.RLIM_INFINITY else min(value, hard)
        resource.setrlimit(which, (min(value, new_hard), new_hard))
    except (ValueError, OSError):
        pass  # best effort; some limits are not settable on every platform


def _drain(stream, cap: int, sink: list[bytes], flags: dict[str, bool]) -> None:
    """Read ``stream`` into ``sink`` up to ``cap`` bytes, then keep draining."""
    total = 0
    while True:
        chunk = stream.read(65536)
        if not chunk:
            break
        if total < cap:
            room = cap - total
            sink.append(chunk[:room])
            total += min(room, len(chunk))
            if total >= cap:
                flags["truncated"] = True
        # past the cap we keep reading (and discarding) so the child does not
        # block on a full pipe.


def _normalize_roots(data_paths: Mapping[str, object] | Iterable[object] | None) -> list[str]:
    if data_paths is None:
        return []
    values = data_paths.values() if isinstance(data_paths, Mapping) else data_paths
    return [str(Path(p).resolve()) for p in values]


def run(
    code: str,
    *,
    data_paths: Mapping[str, object] | Iterable[object] | None = None,
    extra_syspath: Iterable[object] | None = None,
    timeout_s: float = DEFAULT_TIMEOUT_S,
    output_cap: int = DEFAULT_OUTPUT_CAP,
) -> SandboxResult:
    """Execute untrusted ``code`` under isolation and return its captured output.

    ``data_paths`` are the registry-resolved paths the code is allowed to read
    (a mapping ``{name: path}`` or a plain iterable of paths); all other reads
    and every write are refused. ``extra_syspath`` is prepended to the child's
    ``sys.path`` so an allowlisted installed package (e.g. ``pygenogrove``) can be
    imported despite the child running with ``-S``; the import guard still blocks
    every non-allowlisted import, so this does not widen the boundary. See the
    module docstring for the full set of guarantees and the residual risk.
    """
    roots = _normalize_roots(data_paths)
    syspath = _normalize_roots(extra_syspath)
    script = _build_script(code, roots, syspath)

    with tempfile.TemporaryDirectory(prefix="ggask-sandbox-") as tmp:
        script_path = Path(tmp) / "query.py"
        script_path.write_text(script, encoding="utf-8")

        proc = subprocess.Popen(
            [sys.executable, "-I", "-S", str(script_path)],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=tmp,
            env=_child_env(),
            preexec_fn=(lambda: _apply_limits(timeout_s)) if os.name == "posix" else None,
        )

        out_chunks: list[bytes] = []
        err_chunks: list[bytes] = []
        flags = {"truncated": False}
        readers = [
            threading.Thread(target=_drain, args=(proc.stdout, output_cap, out_chunks, flags)),
            threading.Thread(target=_drain, args=(proc.stderr, output_cap, err_chunks, flags)),
        ]
        for t in readers:
            t.start()

        timed_out = False
        try:
            proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            timed_out = True
            _kill(proc)
            proc.wait()
        finally:
            for t in readers:
                t.join()

    stdout = b"".join(out_chunks).decode("utf-8", "replace")
    stderr = b"".join(err_chunks).decode("utf-8", "replace")
    if timed_out:
        note = f"\n[sandbox] killed after exceeding the {timeout_s:g}s wall-clock limit"
        stderr = (stderr + note) if stderr else note.lstrip()
    return SandboxResult(
        stdout=stdout,
        stderr=stderr,
        returncode=proc.returncode if proc.returncode is not None else -1,
        timed_out=timed_out,
        truncated=flags["truncated"],
    )


def _kill(proc: subprocess.Popen) -> None:
    """Kill the child's whole process group (best effort)."""
    try:
        if os.name == "posix":
            os.killpg(os.getpgid(proc.pid), 9)
        else:  # pragma: no cover - non-POSIX
            proc.kill()
    except (ProcessLookupError, PermissionError, OSError):
        try:
            proc.kill()
        except OSError:
            pass


# --------------------------------------------------------------------------- #
# Warm worker — a persistent sandboxed subprocess for interactive/batch use.
#
# run() spawns a fresh interpreter per query, so an expensive first touch (opening
# a large grove: GroveView reads a whole block directory, ~200ms) is paid every
# time. A Worker keeps ONE guarded interpreter alive and caches GroveView.open by
# path, so that cost is paid once per session and each subsequent query reuses the
# open view (~sub-ms). Same in-child guards as run() (import allowlist, scrubbed
# primitives, read-only open, stripped env, memory / no-write rlimits, session
# process group).
#
# Two deliberate differences from run(), because the process is long-lived:
#   * No RLIMIT_CPU — it is cumulative over the process lifetime and would kill the
#     worker after a few queries. The per-query bound is a *wall-clock* deadline
#     enforced by the parent: on overrun the worker is killed and transparently
#     restarted (losing the warm groves). This keeps the hard timeout guarantee.
#   * The worker's own infra modules (io/struct/traceback) are imported before the
#     guard, so they stay importable by query code — pure compute, no os/net/exec,
#     so the boundary is unchanged.
# Protocol: length-prefixed frames over the worker's stdin (request = query code)
# and stdout (response = JSON {stdout, stderr, rc}); query prints are captured into
# buffers, so the real stdout carries only the protocol.
# --------------------------------------------------------------------------- #

# Imported before the guard so the request loop can use them; then GroveView.open is
# wrapped with a per-path cache (the whole point of the warm worker).
_WORKER_INFRA = '''\
import sys as _sys, io as _io, json as _json, struct as _struct, traceback as _tb
_sys.path[:0] = [p for p in _SYSPATH_JSON if p not in _sys.path]
try:
    import pygenogrove as _pg
    _gv_open, _gv_cache = _pg.GroveView.open, {}
    def _cached_gv_open(path, data_offset=0, __o=_gv_open, __c=_gv_cache):
        k = (str(path), data_offset)
        if k not in __c:
            __c[k] = __o(path, data_offset)   # ~200ms once; reused thereafter
        return __c[k]
    _pg.GroveView.open = staticmethod(_cached_gv_open)
except Exception:
    pass
'''

# Runs after the guard/open are installed. Each frame = one query; a fresh globals
# dict per query (no state leak) except the shared GroveView cache.
_WORKER_LOOP = '''\
_IN, _OUT = _sys.stdin.buffer, _sys.__stdout__.buffer

def _read_exact(n):
    b = b""
    while len(b) < n:
        c = _IN.read(n - len(b))
        if not c:
            return None
        b += c
    return b

while True:
    _hdr = _read_exact(4)
    if _hdr is None:
        break
    _payload = _read_exact(_struct.unpack(">I", _hdr)[0])
    if _payload is None:
        break
    _obuf, _ebuf = _io.StringIO(), _io.StringIO()
    _so, _se, _rc = _sys.stdout, _sys.stderr, 0
    _sys.stdout, _sys.stderr = _obuf, _ebuf
    try:
        exec(compile(_payload.decode("utf-8"), "<query>", "exec"), {"__name__": "__main__"})
    except SystemExit:
        pass
    except BaseException:
        _tb.print_exc()
        _rc = 1
    finally:
        _sys.stdout, _sys.stderr = _so, _se
    _resp = _json.dumps({"stdout": _obuf.getvalue()[:_OUTPUT_CAP],
                         "stderr": _ebuf.getvalue()[:_OUTPUT_CAP], "rc": _rc}).encode("utf-8")
    _OUT.write(_struct.pack(">I", len(_resp)) + _resp)
    _OUT.flush()
'''

_MAX_FRAME = 64 * 1024 * 1024  # sanity cap on a response frame length


class _ProtoError(Exception):
    """The worker sent a malformed/truncated frame (or died) — trigger a restart."""


def _build_worker_script(roots: list[str], syspath: list[str], output_cap: int) -> str:
    header = (
        f"_ALLOW_JSON = {json.dumps(sorted(ALLOWED_IMPORTS))}\n"
        f"_ROOTS_JSON = {json.dumps(roots)}\n"
        f"_SYSPATH_JSON = {json.dumps(syspath)}\n"
        f"_OUTPUT_CAP = {int(output_cap)}\n"
    )
    return header + _WORKER_INFRA + _BOOTSTRAP + "\n" + _WORKER_LOOP


def _apply_worker_limits() -> None:  # pragma: no cover - runs in child
    """Pre-exec hook for the long-lived worker: own session + memory / no-write /
    fd caps, but NO RLIMIT_CPU (cumulative — the parent enforces per-query timeouts)."""
    os.setsid()
    if resource is None:
        return
    _setrlimit(resource.RLIMIT_AS, _DEFAULT_MEMORY_BYTES)
    _setrlimit(resource.RLIMIT_FSIZE, 0)
    _setrlimit(resource.RLIMIT_NOFILE, _MAX_OPEN_FDS)


class Worker:
    """A persistent sandboxed interpreter that keeps groves open across queries.

    ``submit(code)`` runs one query and returns a :class:`SandboxResult`, reusing the
    warm groves. On timeout or a dead/broken worker it kills and restarts transparently
    (the next ``submit`` re-opens the groves). ``close()`` tears it down. Same guards as
    :func:`run`; see the section comment for the two long-lived-process differences.
    """

    def __init__(
        self,
        *,
        data_paths: Mapping[str, object] | Iterable[object] | None = None,
        extra_syspath: Iterable[object] | None = None,
        timeout_s: float = DEFAULT_TIMEOUT_S,
        output_cap: int = DEFAULT_OUTPUT_CAP,
    ) -> None:
        self._roots = _normalize_roots(data_paths)
        self._syspath = _normalize_roots(extra_syspath)
        self._timeout_s = timeout_s
        self._output_cap = output_cap
        self._tmp = tempfile.mkdtemp(prefix="ggask-worker-")
        Path(self._tmp, "worker.py").write_text(
            _build_worker_script(self._roots, self._syspath, output_cap), encoding="utf-8"
        )
        self._proc: subprocess.Popen | None = None
        self._start()

    def _start(self) -> None:
        self._proc = subprocess.Popen(
            [sys.executable, "-I", "-S", str(Path(self._tmp, "worker.py"))],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            bufsize=0, cwd=self._tmp, env=_child_env(),
            preexec_fn=_apply_worker_limits if os.name == "posix" else None,
        )

    def _restart(self) -> None:
        if self._proc is not None:
            _kill(self._proc)
            try:
                self._proc.wait(timeout=2)
            except Exception:
                pass
        self._start()

    def submit(self, code: str) -> SandboxResult:
        if self._proc is None or self._proc.poll() is not None:
            self._restart()
        try:
            payload = code.encode("utf-8")
            self._proc.stdin.write(struct.pack(">I", len(payload)) + payload)  # type: ignore[union-attr]
            self._proc.stdin.flush()  # type: ignore[union-attr]
            resp = self._read_frame(self._timeout_s)
        except TimeoutError:
            self._restart()
            return SandboxResult(
                "", f"[worker] killed after exceeding the {self._timeout_s:g}s wall-clock limit",
                -1, timed_out=True)
        except (OSError, _ProtoError, ValueError):
            self._restart()
            return SandboxResult("", "[worker] execution failed; the worker was restarted", -1)
        trunc = len(resp["stdout"]) >= self._output_cap or len(resp["stderr"]) >= self._output_cap
        return SandboxResult(resp["stdout"], resp["stderr"], int(resp["rc"]), truncated=trunc)

    def _read_frame(self, timeout_s: float) -> dict:
        deadline = time.monotonic() + timeout_s
        (n,) = struct.unpack(">I", self._read_n(4, deadline))
        if n > _MAX_FRAME:
            raise _ProtoError("frame too large")
        return json.loads(self._read_n(n, deadline).decode("utf-8"))

    def _read_n(self, n: int, deadline: float) -> bytes:
        fd = self._proc.stdout.fileno()  # type: ignore[union-attr]
        buf = b""
        while len(buf) < n:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError
            if not select.select([fd], [], [], remaining)[0]:
                raise TimeoutError
            chunk = os.read(fd, n - len(buf))
            if not chunk:
                raise _ProtoError("worker closed the pipe")
            buf += chunk
        return buf

    def close(self) -> None:
        if self._proc is not None:
            _kill(self._proc)
            try:
                self._proc.wait(timeout=2)
            except Exception:
                pass
            self._proc = None
        shutil.rmtree(self._tmp, ignore_errors=True)

    def __enter__(self) -> "Worker":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def _selftest() -> None:
    """Protocol/capture/timeout/restart check — stdlib only (no pygenogrove needed)."""
    w = Worker(data_paths=[], extra_syspath=[], timeout_s=2.0)
    try:
        r = w.submit('print("hello"); print("world")')
        assert r.stdout.split() == ["hello", "world"] and r.returncode == 0, r
        r = w.submit('x = 1 / 0')                       # error -> rc 1, traceback in stderr
        assert r.returncode == 1 and "ZeroDivisionError" in r.stderr, r
        r = w.submit('import re; print(re.sub("a", "b", "banana"))')  # warm reuse, allowlisted import
        assert r.stdout.strip() == "bbnbnb", r
        r = w.submit('while True: pass')                # timeout -> kill + restart
        assert r.timed_out, r
        r = w.submit('print("alive after restart")')    # transparently restarted
        assert "alive after restart" in r.stdout, r
        r = w.submit('open("/etc/passwd")')             # guard still holds post-restart
        assert r.returncode == 1 and "restricted" in r.stderr, r
    finally:
        w.close()
    print("worker selftest OK")


if __name__ == "__main__":
    _selftest()
