# SPDX-License-Identifier: GPL-3.0-or-later
"""Host-side enhancer layer: target parsing, gene→TSS resolution, the ENHANCERS injection,
and (when the index bundle is present) the tabix lookups. The generated query code can't read
the index — the sandbox allowlist blocks file I/O — so these run host-side only."""
import json

import pytest

from canopy import llm
from canopy.layers import enhancers

FLAGSHIP = "EFO:0005726"  # LNCaP


def test_parse_cohort_targets_and_code():
    text = ('reasoning...\nCOHORT: MCF-7\n'
            'TARGETS: [{"gene": "MYC"}, {"region": "chr8:127700000-127740000"}]\n'
            "```python\nprint(1)\n```")
    cohort, targets, code = llm.parse_targets_and_code(text)
    assert cohort == "MCF-7"
    assert targets == [{"gene": "MYC"}, {"region": "chr8:127700000-127740000"}]
    assert code == "print(1)\n"


def test_parse_none_declared_is_empty():
    # a structural (non-enhancer) reply declares nothing -> backward compatible
    cohort, targets, code = llm.parse_targets_and_code("```python\nx = 1\n```")
    assert cohort == "" and targets == [] and code == "x = 1\n"


def test_parse_malformed_targets_tolerated():
    cohort, targets, code = llm.parse_targets_and_code("TARGETS: [not json\n```python\np()\n```")
    assert cohort == "" and targets == [] and code == "p()\n"


@pytest.mark.parametrize("text", [
    # declarations in a BARE fence, program in a python fence (the real crash: 'MCF-7' -> NameError)
    '```\nCOHORT: MCF-7\nTARGETS: [{"gene": "EGFR"}]\n```\n```python\np()\n```',
    # declarations leaked INSIDE the python fence as leading lines
    '```python\nCOHORT: MCF-7\nTARGETS: [{"gene": "EGFR"}]\np()\n```',
    # clean: plain declaration lines + one python fence
    'COHORT: MCF-7\nTARGETS: [{"gene": "EGFR"}]\n```python\np()\n```',
])
def test_declarations_never_leak_into_code(text):
    cohort, targets, code = llm.parse_targets_and_code(text)
    assert cohort == "MCF-7" and targets == [{"gene": "EGFR"}]
    assert "COHORT" not in code and "TARGETS" not in code and code.strip() == "p()"


def test_resolve_gene():
    # gene_tss.tsv.gz ships in the package, so this works without the index
    hit = enhancers.resolve_gene("MYC")
    assert hit == ("ENSG00000136997", "chr8", 127735434)
    assert enhancers.resolve_gene("NOT_A_GENE") is None


def test_preamble_roundtrips():
    records = [{"chrom": "chr8", "start": "1", "class": "genic", "score_max": "0.9"}]
    pre = enhancers.preamble(records)
    assert pre.startswith("ENHANCERS = json.loads(")
    loaded = eval(pre.split("=", 1)[1].strip(), {"json": json})  # json is on the sandbox allowlist
    assert loaded == records


@pytest.mark.skipif(not enhancers.index_present(FLAGSHIP),
                    reason="rE2G index bundle not present (download/build it first)")
def test_fetch_for_targets_gene_and_region():
    recs = enhancers.fetch_for_targets(
        [{"gene": "MYC"}, {"region": "chr8:127700000-127740000"}], [FLAGSHIP])
    assert recs
    assert all("target_gene" in r and "score_max" in r for r in recs)
    # deduped across the two overlapping targets
    keys = {(r["chrom"], r["start"], r["end"], r["target_gene"], r["cohort"]) for r in recs}
    assert len(keys) == len(recs)
