# SPDX-License-Identifier: GPL-3.0-or-later
"""Augment the GENCODE grove with an ENCODE-rE2G cohort's enhancer→gene edges.

    python tools/build_re2g_grove.py ENCSR101SLZ ENCSR357SDP ENCSR864JGD   # LNCaP (3 reps)
    python tools/build_re2g_grove.py ENCSR864JGD                            # a single replicate

Resolves the (pinned) GENCODE `.gg`, downloads each replicate's thresholded
element-gene-links BED (once, cached), and writes a *combined* grove: GENCODE
gene/transcript/exon nodes + enhancer nodes + `regulates` edges (score = max over
replicates, n = replicate support) hung on the existing gene keys. Prints the counts
and the final `.gg` size. This is the local build; the produced `.gg` is what gets
pinned to Zenodo and wired into the query path.

Run `--selftest` to exercise the augmenter (multi-replicate merge, edge n, self-promoter,
missed target) on a synthetic 1-gene grove — no network, no htslib, no GENCODE download.
"""

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from genogrove_canopy import resources  # noqa: E402


def _row(**kw):
    """An rE2G edge dict (only the kept columns), overridable by keyword."""
    base = {"chr": "chr8", "start": "0", "end": "0", "class": "intergenic",
            "TargetGene": "GENE", "TargetGeneEnsemblID": "ENSG0", "TargetGeneTSS": "0",
            "isSelfPromoter": "FALSE", "Score": "0.5"}
    base.update(kw)
    return base


def _selftest() -> None:
    """Augment a synthetic grove with TWO cohorts (one of two replicates, one of one) and
    confirm: enhancer + bidirectional regulates/regulated_by edges land on the gene key; a
    link's byCohort map carries per-cohort {score=max, n=replicate support}; a link shared
    across cohorts merges into one edge with both cohorts; gene → its enhancers works via
    regulated_by; an off-target is counted; and it survives a serialize/GroveView round-trip."""
    import pygenogrove as pg

    # A minimal GENCODE-like base grove: one gene node, versioned ENSG id.
    base = pg.Grove(order=100)
    base.insert("chr8", pg.GenomicCoordinate("+", 127735000, 127742000),
                {"type": "gene", "id": "ENSG00000136997.20", "name": "MYC"})

    base_gg = Path(tempfile.mkdtemp()) / "base.gg"
    base.serialize(str(base_gg))

    # Rows target MYC (TSS inside its body) except GHOST (not in the base grove).
    myc = {"TargetGene": "MYC", "TargetGeneEnsemblID": "ENSG00000136997", "TargetGeneTSS": "127736000"}
    cohorts = {
        # "prostate": two replicates. E1 (100-200) in both, E2 (300-400) in rep1 only.
        "prostate": [
            [_row(**myc, start="100", end="200", Score="0.9"),
             _row(**myc, start="300", end="400", Score="0.7"),
             _row(start="500", end="600", TargetGene="GHOST",  # off-target -> missed
                  TargetGeneEnsemblID="ENSG99999999", TargetGeneTSS="127736000")],
            [_row(**myc, start="100", end="200", Score="0.8")],  # E1 again (lower score)
        ],
        # "breast": one replicate. Shares E1 with prostate; adds E3 (700-800).
        "breast": [
            [_row(**myc, start="100", end="200", Score="0.5"),
             _row(**myc, start="700", end="800", Score="0.6")],
        ],
    }
    g, stats = resources.augment_grove(base_gg, cohorts)
    # elements E1,E2,E3 = 3 enhancers; 3 distinct links; GHOST missed once.
    assert stats["enhancers"] == 3 and stats["regulates"] == 3 and stats["missed_targets"] == 1, stats
    # prostate: E1(n=2), E2(n=1) -> 2 links, n_dist {1:1,2:1}; breast: E1(n=1), E3(n=1) -> 2 links.
    assert stats["cohorts"]["prostate"] == {"links": 2, "n_dist": {1: 1, 2: 1}}, stats["cohorts"]
    assert stats["cohorts"]["breast"] == {"links": 2, "n_dist": {1: 2}}, stats["cohorts"]

    out = Path(tempfile.mkdtemp()) / "combined.gg"
    g.serialize(str(out))
    gv = pg.GroveView.open(str(out))

    # E1 (both cohorts): variant ∩ enhancer -> regulates -> MYC; byCohort has both tissues.
    # "*" wildcard strand — enhancers are unstranded, genes are +/-.
    e1 = [h for h in gv.intersect(pg.GenomicCoordinate("*", 150, 150), "chr8")
          if h.data.get("type") == "enhancer"][0]
    assert e1.value.end == 199, e1.value.end  # half-open [100,200) -> closed [100,199]
    fwd = gv.get_edges(e1)[0]
    assert fwd["rel"] == "regulates", fwd
    assert fwd["byCohort"] == {"prostate": {"score": 0.9, "n": 2},   # max(0.9,0.8), 2/2 reps
                               "breast": {"score": 0.5, "n": 1}}, fwd
    targets = gv.get_neighbors_if(e1, lambda m: m.get("rel") == "regulates")
    assert targets[0].data["id"] == "ENSG00000136997.20", targets[0].data  # the GENCODE MYC key

    # Reverse: from the MYC gene node, "its enhancers" via regulated_by = E1, E2, E3.
    myc_key = list(gv.intersect(pg.GenomicCoordinate("*", 127738000, 127738000), "chr8"))
    myc_key = [k for k in myc_key if k.data.get("type") == "gene"][0]
    enhs = gv.get_neighbors_if(myc_key, lambda m: m.get("rel") == "regulated_by")
    assert {e.value.start for e in enhs} == {100, 300, 700}, [e.value.start for e in enhs]
    print("augment_grove selftest OK:", {k: stats[k] for k in ("enhancers", "regulates", "cohorts")})


def main(argv):
    if not argv or argv[0] == "--selftest":
        _selftest()
        return 0
    accessions = sorted(argv)  # a cohort's replicate accessions (one or more)
    label = next((e["biosample_term"] for e in resources.re2g_catalog()
                  if e["accession"] == accessions[0]), "cohort")  # the cohort's name
    print(f"Augmenting GENCODE with rE2G cohort {label!r} {accessions} (GENCODE .gg + edge "
          "download + in-memory deserialize; first run is slow)…", file=sys.stderr)
    base_gg = resources.ensure_all_grove("gencode.human")
    grove, stats = resources.augment_grove(
        base_gg, {label: [resources.re2g_edges(a, "") for a in accessions]})
    out = resources.augmented_grove_path("gencode.human", {label: accessions})  # same cache the CLI uses
    out.parent.mkdir(parents=True, exist_ok=True)
    grove.serialize(str(out))
    print(f"grove: {out}")
    print(f"size:  {out.stat().st_size / 1e6:.1f} MB")
    print(f"stats: {stats['enhancers']} enhancers, {stats['regulates']} regulates edges "
          f"(+ reverse regulated_by), {stats['missed_targets']} rE2G targets not in GENCODE, "
          f"{stats['self_promoters']} self-promoters")
    for name, s in stats["cohorts"].items():
        total = s["links"] or 1
        dist = "  ".join(f"n={k}: {v} ({100 * v / total:.0f}%)" for k, v in s["n_dist"].items())
        print(f"  {name}: {s['links']} links — replicate support: {dist}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
