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
    cohorts, note = cli._resolve_query_cohorts([], "K562; HepG2 ;", ["enhancers"])
    assert seen == [["K562", "HepG2"]] and note is None and list(cohorts) == ["K562", "HepG2"]
    # the user's override wins over the model's term
    cli._resolve_query_cohorts(["LNCaP"], "K562", ["enhancers"])
    assert seen[-1] == ["LNCaP"]


def test_resolve_cohorts_gives_both_layers_keys():
    """One spec -> per-layer keys: a tissue term (bridge row: rE2G id + PCAWG codes), an rE2G
    ontology id or biosample-name substring (rE2G only), a PCAWG project code (SV only)."""
    from genogrove_canopy import cli

    r = cli._resolve_cohorts(["Breast Cancer", "breast", "EFO:0002067", "hepg2", "BRCA-US", "bone"])
    assert r["breast cancer"] == {"re2g": ["EFO:0001203"], "pcawg": ["BRCA-US", "BRCA-UK", "BRCA-EU"], "hichip": ["TCGA-BRCA"]}
    assert r["breast"] == {"re2g": ["UBERON:0008367"], "pcawg": ["BRCA-US", "BRCA-UK", "BRCA-EU"], "hichip": ["TCGA-BRCA"]}
    assert r["K562"] == {"re2g": ["EFO:0002067"], "pcawg": [], "hichip": []}
    assert r["HepG2"]["re2g"] == ["EFO:0001187"] and r["HepG2"]["pcawg"] == []
    assert r["BRCA-US"] == {"re2g": [], "pcawg": ["BRCA-US"], "hichip": []}
    assert r["bone"] == {"re2g": [], "pcawg": ["BOCA-UK", "SARC-US"], "hichip": []}   # rE2G has no bone biosample
    assert cli._resolve_cohorts(["TCGA-BRCA"])["TCGA-BRCA"] == {"re2g": [], "pcawg": [], "hichip": ["TCGA-BRCA"]}
    assert cli._resolve_cohorts(["adrenal"])["adrenal"]["hichip"] == ["TCGA-ACC"]   # the one HiChIP-only term
    assert cli._cohort_ids(r) == ["EFO:0001203", "UBERON:0008367", "EFO:0002067", "EFO:0001187"]
    with pytest.raises(SystemExit, match="no cohort matches"):
        cli._resolve_cohorts(["nonsense-tissue"])


def test_session_attach_wires_the_layer_for_the_resolved_cohort(monkeypatch, tmp_path):
    """`attach` decides what the host prepares: an SV attach for a PCAWG cohort extracts that
    cohort's table, and the next `run`'s script attaches it and names it in SV_COHORTS; no
    enhancer index is touched. Before any attach, `run` uses the base preamble."""
    from genogrove_canopy import cli, preamble
    from genogrove_canopy.layers import enhancers, sv

    table = tmp_path / "BRCA-US.tsv"
    table.write_text("\t".join(sv._FIELDS + ("sample", "cohort")) + "\n")
    monkeypatch.setattr(sv, "cohort_file", lambda code: table)
    monkeypatch.setattr(enhancers, "ensure_index", lambda c: (_ for _ in ()).throw(AssertionError("enhancers not declared")))
    scripts = []

    class _R:
        truncated = False
        returncode, timed_out, stdout, stderr = 0, False, "ok: 1", ""
    base = preamble.build("/x.gg")                    # the real base, as _grove_context hands it over
    session = cli.Session(gg="/x.gg", base=base, execute=lambda s: scripts.append(s) or _R(), say=lambda m: None)
    assert session.run("print(0)\n") == "ok: 1" and "globals().pop('_CANOPY_STATE', None)" in scripts[0]
    note = session.attach(["sv"], "BRCA-US")
    assert 'SV_COHORTS = ["BRCA-US"]' in note and "SVs attached from BRCA-US" in note
    assert session.run("print(1)\n") == "ok: 1" and session.stdout == "ok: 1" and session.code == "print(1)\n"
    script = scripts[1]
    assert "def attach_tracked(" in script and f'["BRCA-US", "{table}"]' in script
    assert 'SV_COHORTS = ["BRCA-US"]' in script and "COHORTS = []" in script
    # exactly one preamble: the base one evicts the worker memo, so it must not precede the
    # cohort one (that concatenation is what silently disabled the warm grove)
    assert script.count("GROVE = ") == 1 and "globals().pop('_CANOPY_STATE', None)" not in script
    assert session.attach([], "") .startswith("error:")
    # the reason nothing resolved reaches the model, not only the user's terminal
    monkeypatch.setattr(cli, "_resolve_cohorts", lambda specs: (_ for _ in ()).throw(SystemExit("no cohort matches")))
    assert "no cohort matched 'nowhere'" in session.attach(["sv"], "nowhere")
    assert "name a tissue" in session.attach(["sv"], "")


def test_session_build_makes_every_later_run_warm(monkeypatch):
    """`build` runs with the mutable grove; after it, a plain `run` must not fall back to the
    lazy base preamble (which would not see the inserted records) — it takes the memo path with
    the read-only view. A failed build still flips the switch: partial inserts are in the memo."""
    from genogrove_canopy import cli

    scripts = []

    class _R:
        truncated = False
        returncode, timed_out, stderr = 0, False, ""

        def __init__(self, out):
            self.stdout = out
    outs = iter(["built: 3", "rows", ""])
    session = cli.Session(gg="/x.gg", base="BASE\n", execute=lambda s: scripts.append(s) or _R(next(outs)), say=lambda m: None)
    assert session.build("GROVE.insert(...)\n") == "built: 3"
    assert "GROVE = _readonly(_grove, ('insert', 'add_edge', 'add_external_key'))" in scripts[0]
    assert session.run("q()\n") == "rows" and session.built
    assert "GROVE = _readonly(_grove)" in scripts[1] and "BASE" not in scripts[1] and "GroveView.open" not in scripts[1]
    assert session.build("x\n").startswith("(built") and "own data built" in session.notes[0]


def test_own_data_files_are_named_for_build_and_granted(monkeypatch, tmp_path):
    """`--data FILE` reaches the model as a path to `open()` in `build`, and the sandbox roots.
    Without any, the block says so; a non-file is refused before anything runs."""
    from genogrove_canopy import cli, resources

    block = cli.resources_block("GROVE", "d", "- l", ["/tmp/mine.bed"])
    assert "`/tmp/mine.bed`" in block and "`build`" in block
    assert "No own files" in cli.resources_block("GROVE", "d", "- l")
    gg = tmp_path / "g.gg"
    gg.touch()
    mine = tmp_path / "mine.bed"
    mine.write_text("chr1\t1\t2\n")
    monkeypatch.setattr(resources, "ensure_all_grove", lambda name: gg)
    block, _pre, paths = cli._grove_context([str(mine)])
    assert str(mine.resolve()) in block and paths[-1] == str(mine.resolve())
    with pytest.raises(SystemExit, match="not a file"):
        cli._grove_context([str(tmp_path)])


def test_answer_renders_the_models_reply_over_the_last_rows(monkeypatch, tmp_path, capsys):
    """Text format: the model's prose, then the rows of its last successful run. A machine
    format keeps stdout to the rows and sends the prose to stderr, so `--format json` still
    parses. An abstain (empty reply at the turn cap) is an error naming the cap; the trace file
    holds every call. A later successful run that prints nothing clears earlier rows."""
    import json
    from genogrove_canopy import cli, llm

    class _R:
        truncated = False
        returncode, timed_out, stderr = 0, False, ""

        def __init__(self, out):
            self.stdout = out
    row = '{"chrom": "chr8", "start": 1, "end": 2, "name": "MYC"}\n'
    outs = iter([row, row, row, row, ""])

    def fake_answer(question, system_prompt, session, **kw):
        session.run("a()\n")
        return "MYC it is.", [{"tool": "run", "input": {"code": "a()\n"}, "output": "x"}]
    monkeypatch.setattr(llm, "answer", fake_answer)
    args = type("A", (), {"model": "m", "show_code": False, "cohort": None, "format": "text",
                          "trace": str(tmp_path / "t.json")})()
    execute = lambda s: _R(next(outs))
    out, err, trace, _, built = cli._answer("q", system_prompt="", base="", gg="/x.gg", args=args, execute=execute)
    assert err == "" and out.startswith("MYC it is.\n\n") and "MYC" in out.splitlines()[-1] and not built
    assert json.loads((tmp_path / "t.json").read_text())["trace"] == trace
    args.trace = str(tmp_path / "no-such-dir" / "t.json")   # an unwritable trace path never costs the answer
    out, err, *_ = cli._answer("q", system_prompt="", base="", gg="/x.gg", args=args, execute=execute)
    assert out.startswith("MYC it is.") and err == "" and "--trace" in capsys.readouterr().err
    args.trace = str(tmp_path / "t.json")
    args.format = "json"
    out, err, *_ = cli._answer("q", system_prompt="", base="", gg="/x.gg", args=args, execute=execute)
    assert json.loads(out)["name"] == "MYC" and "MYC it is." in capsys.readouterr().err

    def two_runs(question, system_prompt, session, **kw):
        assert session.run("a()\n").startswith("{")
        assert session.run("b()\n") == "(the program produced no output)"
        return "Nothing matched.", []
    monkeypatch.setattr(llm, "answer", two_runs)
    out, err, *_ = cli._answer("q", system_prompt="", base="", gg="/x.gg", args=args, execute=execute)
    assert out == "" and err == ""                                      # no stale rows from run a
    monkeypatch.setattr(llm, "answer", lambda *a, **kw: ("", [{"tool": "run"}] * 2))
    out, err, *_ = cli._answer("q", system_prompt="", base="", gg="/x.gg", args=args, execute=None)
    assert out == "" and "2 tool calls" in err and str(llm.MAX_TURNS) in err


def test_sv_question_without_a_cohort_gets_no_default(monkeypatch):
    """The LNCaP enhancer default must not reach an SV question: there is no honest default
    tumour cohort, so resolution returns nothing and says why."""
    from genogrove_canopy import cli

    cohorts, note = cli._resolve_query_cohorts([], "", ["sv"])
    assert cohorts == {} and "name a tissue or pass --cohort" in note
    cohorts, note = cli._resolve_query_cohorts([], "", ["enhancers", "sv"])
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
    cohorts = {"K562": {"re2g": ["EFO:0002067"], "pcawg": [], "hichip": []}}
    cohort_links, sv_files, hichip_files = cli.prepare_layers(cohorts, ["enhancers", "sv", "hichip"], said.append)
    assert cohort_links == {"EFO:0002067": str(links)} and sv_files == {} and hichip_files == {}
    assert any(m.startswith("no PCAWG cohort for cohort(s) K562") for m in said)
    assert any(m.startswith("no TCGA HiChIP cohort for cohort(s) K562") for m in said)


def test_bridge_table_keys_exist_in_both_catalogs():
    """`data/cohorts.tsv` is curated by hand: every rE2G id must be a real biosample, every
    PCAWG code a real cohort, no term or alias may appear twice, and every PCAWG cohort must be
    reachable from some term — otherwise a tumour cohort exists that no tissue word can load."""
    from genogrove_canopy import cli, resources

    rows = cli._bridge()
    re2g = {c["ontology_id"] for c in resources.re2g_cohorts()}
    pcawg = {r["project_code"] for r in resources.pcawg_cohorts()}
    hichip = {r["project_code"] for r in resources.hichip_cohorts()}
    words = [w.lower() for r in rows for w in [r["term"], *r["aliases"]]]
    assert len(words) == len(set(words)), "duplicate term/alias"
    assert all(i in re2g for r in rows for i in r["re2g"])
    assert all(c in pcawg for r in rows for c in r["pcawg"])
    assert pcawg == {c for r in rows for c in r["pcawg"]}, "a PCAWG cohort no term reaches"
    assert all(c in hichip for r in rows for c in r["hichip"])
    assert hichip == {c for r in rows for c in r["hichip"]}, "a HiChIP cohort no term reaches"
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
