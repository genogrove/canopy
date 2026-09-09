# SPDX-License-Identifier: GPL-3.0-or-later
"""Smoke tests for the CLI skeleton."""

import pytest

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


def test_declared_cohorts_split_on_semicolon(monkeypatch):
    """`COHORT: K562; HepG2` — a comparison question names several cohorts on one line."""
    from genogrove_canopy import cli

    seen = []
    monkeypatch.setattr(cli, "_resolve_cohorts", lambda specs: seen.append(specs) or {s: [] for s in specs})
    args = type("A", (), {"cohort": None})()
    cohorts, note = cli._resolve_query_cohorts(args, "K562; HepG2 ;")
    assert seen == [["K562", "HepG2"]] and note is None and list(cohorts) == ["K562", "HepG2"]


def test_resolve_cohorts_gives_both_layers_keys():
    """One spec -> per-layer keys: a tissue term (bridge row: rE2G id + PCAWG codes), an rE2G
    ontology id or biosample-name substring (rE2G only), a PCAWG project code (SV only)."""
    from genogrove_canopy import cli

    r = cli._resolve_cohorts(["Breast Cancer", "EFO:0002067", "hepg2", "BRCA-US", "bone"])
    assert r["breast"] == {"re2g": ["EFO:0001203"], "pcawg": ["BRCA-US", "BRCA-UK", "BRCA-EU"]}
    assert r["K562"] == {"re2g": ["EFO:0002067"], "pcawg": []}
    assert r["HepG2"]["re2g"] == ["EFO:0001187"] and r["HepG2"]["pcawg"] == []
    assert r["BRCA-US"] == {"re2g": [], "pcawg": ["BRCA-US"]}
    assert r["bone"] == {"re2g": [], "pcawg": ["BOCA-UK", "SARC-US"]}   # rE2G has no bone biosample
    assert cli._cohort_ids(r) == ["EFO:0001203", "EFO:0002067", "EFO:0001187"]
    with pytest.raises(SystemExit, match="no cohort matches"):
        cli._resolve_cohorts(["nonsense-tissue"])
