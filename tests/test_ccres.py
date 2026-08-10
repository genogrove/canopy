# SPDX-License-Identifier: GPL-3.0-or-later
"""cCRE node layer: the host-side tabix read and insert-as-grove-nodes attach. The coordinate
conversions (BED half-open -> grove 0-based closed; tabix 1-based) are the risk, so the tests
pin them against a known element — the dELS at chr7:55,191,712-55,191,899 that contains the
EGFR example variant chr7:55,191,822."""
import shutil
import subprocess

import pytest

from canopy.layers import ccres

pg = pytest.importorskip("pygenogrove")


def test_attach_inserts_intersect_queryable_ccre_nodes():
    records = [{"chrom": "chr7", "start": "55191712", "end": "55191899",
                "rdhs": "EH38D5856940", "accession": "EH38E3775532", "ccre_class": "dELS"}]
    g = pg.Grove(order=100)
    assert ccres.attach(g, records) == 1
    # the variant falls in the element -> intersect returns the node with class verbatim. cCREs
    # are identified by `source`, not `type`: the SO term is the generic regulatory_region, so a
    # dELS is never mislabelled as a confirmed `enhancer`.
    hits = [k.data for k in g.intersect(pg.GenomicCoordinate("*", 55191822, 55191822), "chr7")
            if k.data.get("source") == "ENCODE-SCREEN"]
    assert hits == [{"type": "regulatory_region", "source": "ENCODE-SCREEN", "class": "dELS",
                     "id": "EH38E3775532", "rdhs": "EH38D5856940"}]
    # BED end is half-open -> stored as 0-based closed 55191898, so 55191899 is NOT covered
    assert not [k for k in g.intersect(pg.GenomicCoordinate("*", 55191899, 55191899), "chr7")
                if k.data.get("source") == "ENCODE-SCREEN"]


@pytest.mark.skipif(shutil.which("bgzip") is None or shutil.which("tabix") is None,
                    reason="htslib (bgzip/tabix) not installed")
def test_in_region_reads_only_the_locus(tmp_path, monkeypatch):
    bed = tmp_path / "ccre.bed"
    bed.write_text("chr7\t55189954\t55190136\tEH38D1\tEH38E1\tpELS\n"
                   "chr7\t55191712\t55191899\tEH38D2\tEH38E2\tdELS\n"
                   "chr8\t100\t200\tEH38D3\tEH38E3\tPLS\n")
    gz = tmp_path / "ccre.bed.gz"
    subprocess.run(f"sort -k1,1 -k2,2n {bed} | bgzip -c > {gz}", shell=True, check=True)
    subprocess.run(["tabix", "-p", "bed", str(gz)], check=True)
    monkeypatch.setattr(ccres.resources, "indexed_path", lambda name: gz)

    recs = ccres.in_region("chr7", 55191712, 55191899)  # 0-based closed, over the dELS only
    assert [r["accession"] for r in recs] == ["EH38E2"]
    assert recs[0]["ccre_class"] == "dELS"
    assert set(recs[0]) == set(ccres._FIELDS)


def test_all_records_streams_the_registry():
    if not ccres.resources.is_indexed(ccres._NAME):
        pytest.skip("cCRE index not cached")
    r = next(iter(ccres.all_records()))  # generator — just pull the first record
    assert set(r) == set(ccres._FIELDS)
    assert r["chrom"].startswith("chr") and r["ccre_class"]


def test_baked_grove_returns_ccres_alongside_genes():
    """End-to-end: cCREs are baked into the shipped grove and a single GroveView.intersect at the
    EGFR example variant returns EGFR *and* the dELS cCRE. Skipped until the grove is built."""
    from canopy import resources

    gg = resources._baked_grove_gg("gencode.human")
    if not gg.exists():
        pytest.skip("baked grove not built (run `canopy --init`)")
    hits = [k.data for k in pg.GroveView.open(str(gg))
            .intersect(pg.GenomicCoordinate("*", 55191821, 55191821), "chr7")]
    assert "gene" in {d.get("type") for d in hits}
    ccre = next(d for d in hits if d.get("source") == "ENCODE-SCREEN")
    assert ccre["type"] == "regulatory_region" and ccre["class"] == "dELS"
