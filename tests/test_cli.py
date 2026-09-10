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
    cohorts, note = cli._resolve_query_cohorts(args, "K562; HepG2 ;", ["enhancers"])
    assert seen == [["K562", "HepG2"]] and note is None and list(cohorts) == ["K562", "HepG2"]


def test_resolve_cohorts_gives_both_layers_keys():
    """One spec -> per-layer keys: a tissue term (bridge row: rE2G id + PCAWG codes), an rE2G
    ontology id or biosample-name substring (rE2G only), a PCAWG project code (SV only)."""
    from genogrove_canopy import cli

    r = cli._resolve_cohorts(["Breast Cancer", "breast", "EFO:0002067", "hepg2", "BRCA-US", "bone"])
    assert r["breast cancer"] == {"re2g": ["EFO:0001203"], "pcawg": ["BRCA-US", "BRCA-UK", "BRCA-EU"]}
    assert r["breast"] == {"re2g": ["UBERON:0008367"], "pcawg": ["BRCA-US", "BRCA-UK", "BRCA-EU"]}
    assert r["K562"] == {"re2g": ["EFO:0002067"], "pcawg": []}
    assert r["HepG2"]["re2g"] == ["EFO:0001187"] and r["HepG2"]["pcawg"] == []
    assert r["BRCA-US"] == {"re2g": [], "pcawg": ["BRCA-US"]}
    assert r["bone"] == {"re2g": [], "pcawg": ["BOCA-UK", "SARC-US"]}   # rE2G has no bone biosample
    assert cli._cohort_ids(r) == ["EFO:0001203", "UBERON:0008367", "EFO:0002067", "EFO:0001187"]
    with pytest.raises(SystemExit, match="no cohort matches"):
        cli._resolve_cohorts(["nonsense-tissue"])


def test_answer_wires_the_declared_layers_for_the_resolved_cohort(monkeypatch, tmp_path):
    """LAYERS decides what the host prepares: an SV question for a PCAWG cohort extracts that
    cohort's table, and the script handed to the sandbox attaches it and names it in SV_COHORTS;
    no enhancer index is touched."""
    from genogrove_canopy import cli, llm
    from genogrove_canopy.layers import enhancers, sv

    monkeypatch.setattr(llm, "generate_query", lambda q, sp, model=None: ("BRCA-US", ["sv"], "print(1)\n"))
    table = tmp_path / "BRCA-US.tsv"
    table.write_text("\t".join(sv._FIELDS + ("sample", "cohort")) + "\n")
    monkeypatch.setattr(sv, "cohort_file", lambda code: table)
    monkeypatch.setattr(enhancers, "ensure_index", lambda c: (_ for _ in ()).throw(AssertionError("enhancers not declared")))
    scripts = []

    class _R:
        truncated = False
        returncode, timed_out, stdout, stderr = 0, False, "ok: 1", ""
    args = type("A", (), {"model": "m", "show_code": False, "cohort": None, "format": "tsv"})()
    from genogrove_canopy import preamble
    base = preamble.build("/x.gg")                    # the real base, as _grove_context hands it over
    _out, err, *_ = cli._answer("SVs near MYC in BRCA-US?", system_prompt="", base=base, gg="/x.gg",
                               args=args, execute=lambda s: scripts.append(s) or _R())
    assert err == "" and len(scripts) == 1
    script = scripts[0]
    assert "def attach_tracked(" in script and f'["BRCA-US", "{table}"]' in script
    assert 'SV_COHORTS = ["BRCA-US"]' in script and "COHORTS = []" in script
    # exactly one preamble: the base one evicts the worker memo, so it must not precede the
    # cohort one (that concatenation is what silently disabled the warm grove)
    assert script.count("GROVE = ") == 1 and "globals().pop('_CANOPY_STATE', None)" not in script


def test_sv_question_without_a_cohort_gets_no_default(monkeypatch):
    """The LNCaP enhancer default must not reach an SV question: there is no honest default
    tumour cohort, so resolution returns nothing and says why."""
    from genogrove_canopy import cli

    args = type("A", (), {"cohort": None})()
    cohorts, note = cli._resolve_query_cohorts(args, "", ["sv"])
    assert cohorts == {} and "name a tissue or pass --cohort" in note
    cohorts, note = cli._resolve_query_cohorts(args, "", ["enhancers", "sv"])
    assert cohorts and note == "default"                      # an enhancer question still defaults


def test_prepare_layers_reports_a_declared_layer_with_no_data(monkeypatch, tmp_path):
    """`COHORT: K562` / `LAYERS: enhancers; sv`: enhancers attach, but K562 has no PCAWG cohort —
    that must be said, not left as an empty SV_COHORTS and a plausible '0 rearrangements'."""
    from genogrove_canopy import cli
    from genogrove_canopy.layers import enhancers

    links = tmp_path / "k562.tsv"
    links.write_text("")
    monkeypatch.setattr(enhancers, "ensure_index", lambda c: True)
    monkeypatch.setattr(enhancers, "links_file", lambda c: links)
    said = []
    cohorts = {"K562": {"re2g": ["EFO:0002067"], "pcawg": []}}
    cohort_links, sv_files = cli.prepare_layers(cohorts, ["enhancers", "sv"], said.append)
    assert cohort_links == {"EFO:0002067": str(links)} and sv_files == {}
    assert any(m.startswith("no PCAWG cohort for cohort(s) K562") for m in said)


def test_bridge_table_keys_exist_in_both_catalogs():
    """`data/cohorts.tsv` is curated by hand: every rE2G id must be a real biosample, every
    PCAWG code a real cohort, no term or alias may appear twice, and every PCAWG cohort must be
    reachable from some term — otherwise a tumour cohort exists that no tissue word can load."""
    from genogrove_canopy import cli, resources

    rows = cli._bridge()
    re2g = {c["ontology_id"] for c in resources.re2g_cohorts()}
    pcawg = {r["project_code"] for r in resources.pcawg_cohorts()}
    words = [w.lower() for r in rows for w in [r["term"], *r["aliases"]]]
    assert len(words) == len(set(words)), "duplicate term/alias"
    assert all(i in re2g for r in rows for i in r["re2g"])
    assert all(c in pcawg for r in rows for c in r["pcawg"])
    assert pcawg == {c for r in rows for c in r["pcawg"]}, "a PCAWG cohort no term reaches"
    # the plain tissue word means the tissue biosample where rE2G has one; the disease word the line
    assert cli._resolve_cohorts(["liver"])["liver"]["re2g"] == ["UBERON:0002107"]
    assert cli._resolve_cohorts(["HCC"])["liver cancer"]["re2g"] == ["EFO:0001187"]


@pytest.mark.parametrize("spec", ["", " ", "\t"])
def test_empty_explicit_cohort_is_rejected(spec):
    from genogrove_canopy.cli import _resolve_cohorts
    with pytest.raises(SystemExit, match="empty cohort"):
        _resolve_cohorts([spec])


def test_every_catalogue_name_resolves_to_its_own_biosample():
    from genogrove_canopy import cli, resources
    for cohort in resources.re2g_cohorts():
        result = cli._resolve_cohorts([" " + cohort["name"].swapcase() + " "])
        assert cli._cohort_ids(result) == [cohort["ontology_id"]], cohort["name"]


def test_bridge_words_that_are_biosample_names_name_that_biosample_and_keep_svs():
    """Exact biosample names win over the bridge. So a bridge term or alias that is ALSO an rE2G
    biosample name must point at that very biosample — otherwise the word resolves to the
    biosample, the bridge row is bypassed, and its PCAWG cohorts are silently lost (this was
    `lung`, `pancreas` and `RCC`)."""
    from genogrove_canopy import cli, resources

    names = {c["name"].strip().lower(): c["ontology_id"] for c in resources.re2g_cohorts()}
    for r in cli._bridge():
        for w in [r["term"], *r["aliases"]]:
            if w.lower() in names:
                assert r["re2g"] == [names[w.lower()]], f"{w!r}: bridge row names another biosample"
                assert cli._pcawg_codes(cli._resolve_cohorts([w])) == r["pcawg"], w
