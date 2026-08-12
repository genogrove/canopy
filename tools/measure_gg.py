# SPDX-License-Identifier: GPL-3.0-or-later
"""Measure the whole-genome GENCODE .gg: build time + serialized size vs the GFF.

The number that decides the substrate question — is a prebuilt .gg a viable single
download (serve everything via GroveView), or is it so much bigger than the tabix GFF
that region-tabix stays the primary path?

    uv run python tools/measure_gg.py

Builds via resources.ensure_all_grove (reads the whole cached GFF, region="", serializes
+ caches the _all.gg). A couple of minutes + notable memory the first time; instant after.
"""

from __future__ import annotations

import os
import time

from genogrove_canopy import resources as r

mb = lambda p: os.path.getsize(p) / 1e6

gff = r.indexed_path("gencode.human")               # cached bgzip+tabix GFF
tbi = str(gff) + ".tbi"

t = time.time()
gg = r.ensure_all_grove("gencode.human")            # build whole-genome grove -> .gg (cached)
dt = time.time() - t

print(f"GFF (.gff.gz) : {mb(gff):8.1f} MB   .tbi: {os.path.getsize(tbi)/1e3:6.0f} KB")
print(f".gg (grove)   : {mb(gg):8.1f} MB   built in {dt:5.0f}s")
print(f"ratio .gg/GFF : {mb(gg) / mb(gff):.2f}x")
print(f".gg path      : {gg}")

# Bonus: block count — how many blocks a GroveView pages over (lazy-read granularity).
try:
    import pygenogrove as pg
    v = pg.GroveView.open(str(gg))
    print(f"blocks        : {v.block_count()} (GroveView pages only the ones a query touches)")
except Exception as e:  # pragma: no cover
    print(f"(block_count skipped: {e})")
