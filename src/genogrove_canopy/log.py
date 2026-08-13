# SPDX-License-Identifier: GPL-3.0-or-later
"""Timestamped progress messages for the CLI.

Everything goes to **stderr**, so `canopy … > out.tsv` keeps the answer clean and a pipeline
never has to parse around progress chatter.

Deliberately not `logging`: there is no configuration to expose, no levels anyone would set, and
no library consumer — the CLI wants one line shape and a spinner-free progress counter.
"""

from __future__ import annotations

import sys
import time

_PREFIX = "CANOPY"


def _stamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def say(message: str) -> None:
    """Print one timestamped line: ``[CANOPY - 2026-08-13 09:41:02] message``."""
    print(f"[{_PREFIX} - {_stamp()}] {message}", file=sys.stderr, flush=True)


def took(message: str, seconds: float) -> None:
    """Print a line with an elapsed time, e.g. ``… done (6.9s)``."""
    say(f"{message} ({_fmt_duration(seconds)})")


def _fmt_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    m, s = divmod(int(seconds), 60)
    return f"{m}m{s:02d}s"


def _fmt_size(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if abs(n) < 1024 or unit == "GB":
            return f"{n:.0f}{unit}" if unit in ("B", "KB") else f"{n:.1f}{unit}"
        n /= 1024
    raise AssertionError("unreachable: the GB branch always returns")


class Progress:
    """Percentage counter for a download, redrawn in place on a TTY.

    On a non-TTY (a pipe, CI logs, `2> file`) in-place redraws would produce thousands of junk
    lines, so there it prints one line per decile instead. Either way the final line records the
    size and elapsed time, so a log kept from a non-interactive run still says what happened.
    """

    def __init__(self, label: str, total: int | None) -> None:
        self.label = label
        self.total = total or 0
        self.done = 0
        self.start = time.perf_counter()
        self._tty = sys.stderr.isatty()
        self._last_draw = 0.0
        self._last_decile = -1

    def advance(self, n: int) -> None:
        self.done += n
        now = time.perf_counter()
        if self._tty:
            if now - self._last_draw < 0.1:  # ~10 fps is plenty; redrawing per 1 MB chunk flickers
                return
            self._last_draw = now
            print(f"\r[{_PREFIX} - {_stamp()}] {self.label} {self._bar()}",
                  end="", file=sys.stderr, flush=True)
        elif self.total:
            decile = int(self.done * 10 / self.total)
            if decile > self._last_decile:
                self._last_decile = decile
                say(f"{self.label} {self._bar()}")

    def _bar(self) -> str:
        if not self.total:  # server sent no Content-Length
            return f"{_fmt_size(self.done)}"
        pct = min(100, self.done * 100 // self.total)
        return f"{pct:3d}%  {_fmt_size(self.done)} / {_fmt_size(self.total)}"

    def _clear(self) -> None:
        """End the in-place line. \033[K erases to end of line — without it a shorter following
        line leaves the tail of the longer progress line behind."""
        if self._tty:
            print("\r\033[K", end="", file=sys.stderr, flush=True)

    def finish(self) -> None:
        elapsed = time.perf_counter() - self.start
        rate = self.done / elapsed if elapsed else 0
        self._clear()
        took(f"{self.label} done — {_fmt_size(self.done)} at {_fmt_size(rate)}/s", elapsed)

    def abort(self) -> None:
        """Close out an interrupted transfer.

        The in-place line ends without a newline by design, so whatever prints next lands *on* it.
        Without this, a dropped connection writes its error across the half-drawn progress bar —
        mangling the message exactly when the user most needs to read it.
        """
        self._clear()
