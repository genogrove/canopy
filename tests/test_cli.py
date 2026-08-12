# SPDX-License-Identifier: GPL-3.0-or-later
"""Smoke tests for the CLI skeleton."""

from genogrove_canopy.cli import build_parser, main


def test_parser_builds():
    parser = build_parser()
    args = parser.parse_args(["a question", "--show-code"])
    assert args.question == "a question"
    assert args.show_code is True
    assert args.model == "claude-opus-4-8"


def test_no_question_prints_help_and_succeeds(capsys):
    assert main([]) == 0
    out = capsys.readouterr().out
    assert "canopy" in out
