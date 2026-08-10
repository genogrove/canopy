# SPDX-License-Identifier: GPL-3.0-or-later
"""GFF -> universal Grove load path. Runs only where pygenogrove is installed
(CI); skipped in the bare skeleton env."""

from __future__ import annotations

import pytest

pg = pytest.importorskip("pygenogrove")

from canopy.gff import load_gff

GFF3 = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG1\n"
    "chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tID=exon1;Parent=ENSG1\n"
    "chr2\tHAVANA\tgene\t5000\t6000\t.\t-\t.\tID=ENSG2\n"
)


def _write(tmp_path):
    p = tmp_path / "mini.gff3"
    p.write_text(GFF3)
    return p


def test_type_filter_keeps_only_genes(tmp_path) -> None:
    p = _write(tmp_path)
    assert load_gff(p, types={"gene"}).size() == 2  # the exon is filtered out
    assert load_gff(p).size() == 3


def test_intersect_finds_overlapping_gene(tmp_path) -> None:
    genes = load_gff(_write(tmp_path), types={"gene"})
    inside = pg.GenomicCoordinate("*", 1400, 1600)  # within the chr1 gene
    hits = list(genes.intersect(inside, "chr1"))
    assert len(hits) == 1
    assert hits[0].data["type"] == "gene"  # JSON payload on the universal Grove
    assert hits[0].data["id"] == "ENSG1"
    outside = pg.GenomicCoordinate("*", 3000, 3100)  # past the chr1 gene
    assert len(list(genes.intersect(outside, "chr1"))) == 0


def test_seqid_filter(tmp_path) -> None:
    genes = load_gff(_write(tmp_path), types={"gene"}, seqids={"chr2"})
    assert genes.size() == 1


def test_region_filter(tmp_path) -> None:
    p = _write(tmp_path)  # chr1 gene 1-based [1000,2000] -> 0-based closed [999,1999]
    assert load_gff(p, region=("chr1", 1400, 1600)).size() == 1  # overlaps the gene
    assert load_gff(p, region=("chr1", 2500, 3000)).size() == 0  # window past every chr1 feature
    assert load_gff(p, region=("chr2", 999, 1999)).size() == 0   # right window, wrong chromosome


HIER = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t1000\t3000\t.\t+\t.\tID=g1;gene_name=AAA;gene_type=protein_coding\n"
    "chr1\tHAVANA\ttranscript\t1000\t3000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tHAVANA\texon\t2900\t3000\t.\t+\t.\tID=e2;Parent=t1\n"
)


def test_hierarchy_edges(tmp_path) -> None:
    p = tmp_path / "hier.gff3"
    p.write_text(HIER)
    g = load_gff(p)  # all feature types, so the hierarchy is intact

    gene = next(
        k for k in g.intersect(pg.GenomicCoordinate("*", 1500, 1500), "chr1")
        if k.data["type"] == "gene"
    )
    assert gene.data["name"] == "AAA"
    assert gene.data["biotype"] == "protein_coding"

    txs = list(g.get_neighbors(gene))  # gene -> transcript
    assert [t.data["type"] for t in txs] == ["transcript"]
    assert g.get_edges(gene)[0] == {"rel": "contains"}  # labelled containment edge

    # transcript -> first (5') exon only; the rest hang off the splice chain. Exons are keyed by
    # position, not id — one key per physical exon per gene — so identify them by coordinates.
    entry = list(g.get_neighbors(txs[0]))
    assert [(x.value.start, x.value.end) for x in entry] == [(999, 1099)]   # e1
    assert g.get_edges(txs[0])[0] == {"rel": "first_exon", "tx": "t1"}

    # splice chain: '+' strand, ascending order e1 -> e2. `tx` names the isoform the edge
    # belongs to — an exon shared by several isoforms has one outgoing `next` per isoform.
    assert [(n.value.start, n.value.end) for n in g.get_neighbors(entry[0])] == [(2899, 2999)]
    assert g.get_edges(entry[0])[0] == {"rel": "next", "tx": "t1"}

    # no CDS in this fixture -> the transcript has no span, so the derived exon range is None
    from canopy.gff import exon_cds

    assert txs[0].data["cds_start"] is None
    assert exon_cds(entry[0], txs[0]) is None
    # no id / biotype / cds on an exon. `name` is None because these exon lines carry no
    # gene_name — real GENCODE repeats it on every line (see SHARED / the dedup test).
    assert entry[0].data == {"type": "exon", "name": None}


CODING = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t1000\t3000\t.\t+\t.\tID=g1;gene_name=AAA;gene_type=protein_coding\n"
    "chr1\tHAVANA\ttranscript\t1000\t3000\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1\n"
    "chr1\tHAVANA\texon\t2800\t3000\t.\t+\t.\tID=e2;Parent=t1\n"
    "chr1\tHAVANA\tCDS\t1100\t1200\t.\t+\t0\tID=c1;Parent=t1\n"
    "chr1\tHAVANA\tCDS\t2800\t2900\t.\t+\t2\tID=c2;Parent=t1\n"
    "chr1\tHAVANA\tfive_prime_UTR\t1000\t1099\t.\t+\t.\tID=u1;Parent=t1\n"
)


def test_cds_folded_into_exons(tmp_path) -> None:
    p = tmp_path / "coding.gff3"
    p.write_text(CODING)
    g = load_gff(p)

    everything = list(g.intersect(pg.GenomicCoordinate("*", 0, 5000), "chr1"))
    # CDS and UTR are folded/derived, never inserted as keys
    assert {k.data["type"] for k in everything} == {"gene", "transcript", "exon"}

    tx = next(k for k in everything if k.data["type"] == "transcript")
    assert (tx.data["cds_start"], tx.data["cds_end"]) == (1099, 2899)  # CDS span, 0-based closed

    # The coding sub-range is per (exon, transcript) — an exon key is shared by every isoform
    # that uses it, coding in one and UTR in another — so it is derived, never stored.
    from canopy.gff import exon_cds

    exons = {k.value.start: k for k in everything if k.data["type"] == "exon"}
    assert "cds" not in exons[999].data                       # not on the exon payload
    assert exon_cds(exons[999], tx) == [1099, 1199]           # 999..1098 is 5' UTR (derived)
    assert exon_cds(exons[2799], tx) == [2799, 2899]          # 2900..2999 is 3' UTR (derived)


MINUS = (
    "##gff-version 3\n"
    "chr3\tHAVANA\tgene\t100\t400\t.\t-\t.\tID=mg\n"
    "chr3\tHAVANA\ttranscript\t100\t400\t.\t-\t.\tID=mt;Parent=mg\n"
    "chr3\tHAVANA\texon\t100\t200\t.\t-\t.\tID=lo;Parent=mt\n"
    "chr3\tHAVANA\texon\t300\t400\t.\t-\t.\tID=hi;Parent=mt\n"
)


SHARDED = (
    "##gff-version 3\n"
    "chr1\tH\tgene\t1000\t2000\t.\t+\t.\tID=g1\n"
    "chr2\tH\tgene\t3000\t4000\t.\t-\t.\tID=g2\n"
)


def test_write_sharded_groves_splits_by_chromosome(tmp_path) -> None:
    from canopy.gff import write_sharded_groves

    p = tmp_path / "two.gff3"
    p.write_text(SHARDED)
    out = tmp_path / "out"
    out.mkdir()

    seqids = write_sharded_groves(p, out)
    assert set(seqids) == {"chr1", "chr2"}
    assert pg.Grove.deserialize(str(out / "_all.gg")).size() == 2   # whole genome
    assert pg.Grove.deserialize(str(out / "chr1.gg")).size() == 1   # one shard, one gene
    assert pg.Grove.deserialize(str(out / "chr2.gg")).size() == 1


def test_splice_chain_is_strand_aware(tmp_path) -> None:
    p = tmp_path / "minus.gff3"
    p.write_text(MINUS)
    g = load_gff(p)
    # 5'->3' on the '-' strand runs high coordinate -> low: hi (299..399) -> lo (99..199).
    hi = next(
        k for k in g.intersect(pg.GenomicCoordinate("*", 349, 349), "chr3")
        if k.data["type"] == "exon"
    )
    assert [(n.value.start, n.value.end) for n in g.get_neighbors(hi)] == [(99, 199)]


# A unified GFF: the GENCODE hierarchy plus a foreign layer (ENCODE-SCREEN cCREs) whose
# column-9 vocabulary is its own. The cCRE classification lives only in those attributes.
UNIFIED = (
    "##gff-version 3\n"
    "chr1\tENCODE-SCREEN\tregulatory_region\t900\t1100\t.\t.\t.\t"
    "ID=EH38E1;class=dELS;rdhs=EH38D1\n"
    "chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=g1;gene_name=G1;gene_type=lncRNA\n"
    "chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=t1;Parent=g1;gene_name=G1\n"
    "chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=e1;Parent=t1;gene_name=G1\n"
)


def test_foreign_layer_keeps_its_attributes(tmp_path) -> None:
    """A non-gene/transcript/exon feature keeps its column-9 attributes + source, so a cCRE
    in a unified GFF survives the load with its `class`. Modelled types keep the lean payload
    (GENCODE's tag/level/havana_* must NOT be dragged along). Both loaders must agree."""
    from canopy.gff import build_grove, load_gff

    p = tmp_path / "unified.gff3"
    p.write_text(UNIFIED)
    want = {"type": "regulatory_region", "source": "ENCODE-SCREEN",
            "id": "EH38E1", "class": "dELS", "rdhs": "EH38D1"}

    for g in (load_gff(p), build_grove(str(p))):
        by_type = {k.data["type"]: k.data
                   for k in g.intersect(pg.GenomicCoordinate("*", 1000, 1100), "chr1")}
        assert by_type["regulatory_region"] == want          # class survives
        assert by_type["gene"] == {"type": "gene", "id": "g1",  # modelled: still lean
                                   "name": "G1", "biotype": "lncRNA"}


# A protein-coding gene with an NMD isoform — GENCODE repeats gene_type on every child line,
# so this is the case where a transcript's own biotype and its gene's genuinely disagree.
NMD = (
    "##gff-version 3\n"
    "chr1\tHAVANA\tgene\t1000\t3000\t.\t+\t.\tID=g1;gene_name=AAA;gene_type=protein_coding\n"
    "chr1\tHAVANA\ttranscript\t1000\t3000\t.\t+\t.\tID=t1;Parent=g1;gene_name=AAA;"
    "gene_type=protein_coding;transcript_type=nonsense_mediated_decay\n"
    "chr1\tHAVANA\texon\t1000\t1100\t.\t+\t.\tID=e1;Parent=t1;gene_name=AAA;"
    "gene_type=protein_coding;transcript_type=nonsense_mediated_decay\n"
)


def test_biotype_is_the_features_own(tmp_path) -> None:
    """A transcript stores `transcript_type`, not its gene's — otherwise an NMD isoform of a
    protein-coding gene reads back as protein_coding and the distinction is unrecoverable. The
    gene's biotype is never lost: it's on the gene node. Exons carry no biotype at all."""
    from canopy.gff import build_grove

    p = tmp_path / "nmd.gff3"
    p.write_text(NMD)

    for g in (load_gff(p), build_grove(str(p))):
        by_type = {k.data["type"]: k.data
                   for k in g.intersect(pg.GenomicCoordinate("*", 1000, 1100), "chr1")}
        assert by_type["transcript"]["biotype"] == "nonsense_mediated_decay"
        assert by_type["gene"]["biotype"] == "protein_coding"   # still there, one edge up
        assert "biotype" not in by_type["exon"]                 # exons have none of their own
        assert by_type["exon"]["name"] == "AAA"                 # but still name the gene


# Two isoforms of one gene sharing a middle exon, plus a DIFFERENT gene whose exon happens to
# have the identical interval. GFF emits one exon line per transcript, so the shared exon appears
# three times — but only two of them are the same exon.
SHARED = (
    "##gff-version 3\n"
    "chr1\tH\tgene\t1000\t4000\t.\t+\t.\tID=g1;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\ttranscript\t1000\t4000\t.\t+\t.\tID=t1;Parent=g1;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\ttranscript\t2000\t4000\t.\t+\t.\tID=t2;Parent=g1;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\texon\t1000\t1100\t.\t+\t.\tID=exon:t1:1;Parent=t1;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\texon\t2000\t2100\t.\t+\t.\tID=exon:t1:2;Parent=t1;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\texon\t2000\t2100\t.\t+\t.\tID=exon:t2:1;Parent=t2;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\texon\t3000\t3100\t.\t+\t.\tID=exon:t1:3;Parent=t1;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\texon\t3000\t3100\t.\t+\t.\tID=exon:t2:2;Parent=t2;gene_id=G1;gene_name=AAA\n"
    "chr1\tH\tgene\t2000\t2100\t.\t+\t.\tID=g2;gene_id=G2;gene_name=BBB\n"
    "chr1\tH\ttranscript\t2000\t2100\t.\t+\t.\tID=t3;Parent=g2;gene_id=G2;gene_name=BBB\n"
    "chr1\tH\texon\t2000\t2100\t.\t+\t.\tID=exon:t3:1;Parent=t3;gene_id=G2;gene_name=BBB\n"
)


def _walk(g, tx):
    """An isoform's exons 5'->3', following only the edges tagged with THIS transcript."""
    tid = tx.data["id"]
    entry = list(g.get_neighbors_if(tx, lambda m: m and m["rel"] == "first_exon" and m["tx"] == tid))
    out, cur = [], (entry[0] if entry else None)
    while cur is not None:
        out.append(cur)
        nxt = list(g.get_neighbors_if(cur, lambda m: m and m["rel"] == "next" and m["tx"] == tid))
        cur = nxt[0] if nxt else None
    return out


def test_exons_dedup_per_gene_and_chains_stay_separate(tmp_path) -> None:
    """GFF repeats an exon line per transcript; those are the same exon and collapse to one key.
    An identical interval in a DIFFERENT gene is a different exon and must not merge. The `tx`
    tag on the chain edges is what keeps each isoform's walk on its own path through a shared key."""
    from canopy.gff import build_grove

    p = tmp_path / "shared.gff3"
    p.write_text(SHARED)

    for g in (load_gff(p), build_grove(str(p))):
        at2000 = [k for k in g.intersect(pg.GenomicCoordinate("*", 2050, 2050), "chr1")
                  if k.data["type"] == "exon"]
        assert len(at2000) == 2                                    # 3 GFF lines -> 2 keys
        assert {k.data["name"] for k in at2000} == {"AAA", "BBB"}   # merged within, not across

        txs = {k.data["id"]: k for k in g.intersect(pg.GenomicCoordinate("*", 2050, 2050), "chr1")
               if k.data["type"] == "transcript"}
        # t1 has 3 exons, t2 has 2 — and both run through the SAME shared key at 1999..2099
        assert [(e.value.start, e.value.end) for e in _walk(g, txs["t1"])] == [
            (999, 1099), (1999, 2099), (2999, 3099)]
        assert [(e.value.start, e.value.end) for e in _walk(g, txs["t2"])] == [
            (1999, 2099), (2999, 3099)]
