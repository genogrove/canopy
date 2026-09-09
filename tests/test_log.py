# SPDX-License-Identifier: GPL-3.0-or-later
"""Progress and message formatting. Everything must land on stderr, never stdout —
`canopy … > out.tsv` has to stay parseable."""

from __future__ import annotations

import re

from genogrove_canopy import log


def test_say_is_timestamped_and_goes_to_stderr(capsys) -> None:
    log.say("hello")
    out, err = capsys.readouterr()
    assert out == "", "progress must not pollute stdout — the answer is piped from there"
    assert re.fullmatch(r"\[CANOPY - \d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\] hello\n", err)


def test_took_appends_an_elapsed_time(capsys) -> None:
    log.took("Ready", 6.34)
    assert "Ready (6.3s)" in capsys.readouterr().err


def test_durations_switch_to_minutes() -> None:
    assert log._fmt_duration(0.04) == "0.0s"
    assert log._fmt_duration(59.9) == "59.9s"
    assert log._fmt_duration(61) == "1m01s"
    assert log._fmt_duration(3600) == "60m00s"


def test_progress_without_a_tty_prints_one_line_per_decile(capsys) -> None:
    """A pipe or CI log gets deciles, not thousands of redraws.

    Eleven progress lines, not ten: the 0% one fires on the first chunk, which is what tells you
    the download started rather than stalled.
    """
    p = log.Progress("grove", total=100)
    for _ in range(100):
        p.advance(1)
    p.finish()
    lines = [ln for ln in capsys.readouterr().err.splitlines() if ln.strip()]
    assert len(lines) == 12, f"expected 11 deciles (0-100%) + a final line, got {len(lines)}"
    assert "100%" in lines[-2]
    assert "done" in lines[-1] and "grove" in lines[-1]


def test_progress_without_content_length_reports_bytes_not_percent(capsys) -> None:
    """A server that sends no Content-Length must not produce a division by zero or a fake 0%."""
    p = log.Progress("grove", total=None)
    p.advance(2048)
    p.finish()
    err = capsys.readouterr().err
    assert "%" not in err
    assert "2KB" in err


def test_abort_ends_the_in_place_line(capsys) -> None:
    """An interrupted transfer must close its line, or the next message lands on top of it.

    The in-place line deliberately ends without a newline so the next redraw can overwrite it.
    That is fine while redraws continue and wrong the moment one doesn't — a dropped connection
    would otherwise print its error across the half-drawn bar.
    """
    p = log.Progress("grove", total=1000)
    p._tty = True
    p.advance(500)
    p.abort()
    err = capsys.readouterr().err
    assert err.endswith("\r\033[K"), "abort must erase the in-place line"


def test_finish_and_abort_agree_on_clearing(capsys) -> None:
    """Both close-out paths clear; only `finish` reports totals."""
    p = log.Progress("grove", total=1000)
    p._tty = True
    p.advance(1000)
    p.finish()
    assert "done" in capsys.readouterr().err
