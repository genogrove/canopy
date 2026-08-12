# SPDX-License-Identifier: GPL-3.0-or-later
"""Fetch the ENCODE-rE2G biosample catalog → a TSV that ships with the package.

The catalog is metadata only: one row per prediction annotation (biosample → accession),
so `genogrove_canopy.resources` can show the agent every biosample and fetch the edge BED for only the
one(s) a query names (lazy axis 1; region tabix is axis 2). Regenerate when ENCODE adds
biosamples.

    python tools/fetch_re2g_catalog.py            # writes src/genogrove_canopy/data/encode_re2g_catalog.tsv

Stdlib only (urllib + json + csv) — no third-party deps, runs in any interpreter.
"""

from __future__ import annotations

import csv
import json
import sys
import urllib.request
from pathlib import Path

ENCODE = "https://www.encodeproject.org"
OUT = Path(__file__).resolve().parent.parent / "src" / "ask" / "data" / "encode_re2g_catalog.tsv"

# One search over all released ENCODE-rE2G prediction annotations, GRCh38. `field=` trims
# the payload to what the catalog needs; limit=all returns every biosample in one call.
SEARCH = (
    f"{ENCODE}/search/?type=Annotation&searchTerm=rE2G&assembly=GRCh38&status=released"
    "&limit=all&format=json"
    "&field=accession&field=annotation_type&field=description"
    "&field=biosample_ontology.term_name&field=biosample_ontology.term_id"
    "&field=biosample_ontology.classification"
)

COLUMNS = ("accession", "biosample_term", "biosample_id", "biosample_type", "model")


def fetch() -> list[dict[str, str]]:
    req = urllib.request.Request(SEARCH, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req) as resp:  # noqa: S310 — fixed ENCODE host
        graph = json.load(resp)["@graph"]
    rows = []
    for a in graph:
        # Keep only the enhancer-gene prediction annotations (skip any stray matches).
        if a.get("annotation_type") != "element gene regulatory interaction predictions":
            continue
        bs = a.get("biosample_ontology") or {}
        desc = (a.get("description") or "").lower()
        rows.append({
            "accession": a["accession"],
            "biosample_term": bs.get("term_name", ""),
            "biosample_id": bs.get("term_id", ""),
            "biosample_type": bs.get("classification", ""),
            # ENCODE-rE2G ships DNase-only everywhere + an Extended model for a few tiers.
            "model": "Extended" if "extended" in desc else "DNase-only",
        })
    rows.sort(key=lambda r: (r["biosample_term"], r["accession"]))
    return rows


def main() -> None:
    rows = fetch()
    if not rows:
        sys.exit("no ENCODE-rE2G annotations returned — check the ENCODE portal / query")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    # Sanity: the prostate biosamples we build the hero on must be present.
    prostate = [r for r in rows if "prostate" in r["biosample_term"].lower()
                or "LNCaP" in r["biosample_term"] or r["biosample_term"] in {"PC-3", "RWPE1", "RWPE2"}]
    print(f"wrote {len(rows)} biosamples to {OUT} ({len(prostate)} prostate)", file=sys.stderr)
    assert any("LNCaP" in r["biosample_term"] for r in rows), "LNCaP missing from catalog"


if __name__ == "__main__":
    main()
