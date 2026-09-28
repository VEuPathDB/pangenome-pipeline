"""05_Cross-genome_annotation_integration - Liftoff transfer filtering and overlap resolution."""

import pytest
from conftest import gff_line

SCRIPTS = "05_Cross-genome_annotation_integration/scripts/"


@pytest.fixture(scope="module")
def good(load_script):
    return load_script(SCRIPTS + "get_good_transfer_ids.py")


@pytest.fixture(scope="module")
def ffa(load_script):
    return load_script(SCRIPTS + "filter_final_annotation.py")


def _gene(gid, start, end, segs, seqname="chr1", source="good_X", **kw):
    gene = {
        "id": gid,
        "seqname": seqname,
        "source": source,
        "start": start,
        "end": end,
        "cds_segments": segs,
        "cds_bp": sum(e - s + 1 for s, e in segs),
        "has_domain": False,
        "description": "",
    }
    gene.update(kw)
    return gene


#
# get_good_transfer_ids.py
#


def test_parse_attributes_is_anchored(good):
    attrs = good.parse_attributes("ref_sequence_ID=0.1;ID=g1;sequence_ID=0.95;coverage=1e-05")
    assert attrs["ID"] == "g1"
    assert attrs["sequence_ID"] == "0.95"
    assert good.to_float(attrs["coverage"]) == pytest.approx(1e-05)


@pytest.mark.parametrize(
    ("attrs", "keep_extra", "expected"),
    [
        ("ID=g;valid_ORF=True;coverage=0.95;sequence_ID=0.95", False, True),
        ("ID=g;valid_ORFs=1;coverage=0.95;sequence_ID=0.95", False, True),
        ("ID=g;valid_ORF=False;coverage=1;sequence_ID=1", False, False),
        ("ID=g;coverage=1;sequence_ID=1", False, False),  # no ORF info -> not good
        ("ID=g;valid_ORF=True;coverage=0.85;sequence_ID=1", False, False),
        ("ID=g;valid_ORF=True;coverage=1;sequence_ID=1;extra_copy_number=1", False, False),
        ("ID=g;valid_ORF=True;coverage=1;sequence_ID=1;extra_copy_number=1", True, True),
    ],
)
def test_is_good_transfer(good, attrs, keep_extra, expected):
    assert good.is_good_transfer(good.parse_attributes(attrs), 0.9, 0.9, keep_extra) is expected


def test_extract_good_transfer_ids_only_counts_gene_features(good, write_file):
    text = (
        gff_line("chr1", "protein_coding_gene", 1, 100, "+", "ID=g1;valid_ORF=True;coverage=1;sequence_ID=1")
        + gff_line("chr1", "mRNA", 1, 100, "+", "ID=m1;Parent=g1;valid_ORF=True;coverage=1;sequence_ID=1")
        + gff_line("chr1", "gene", 200, 300, "+", "ID=g2;valid_ORF=True;coverage=0.5;sequence_ID=1")
    )
    ids, n = good.extract_good_transfer_ids(
        write_file("lift.gff3", text), 0.9, 0.9, False, {"gene", "protein_coding_gene"}
    )
    assert ids == ["g1"]
    assert n == 2


#
# filter_final_annotation.py
#


def test_segments_overlap_bp(ffa):
    assert ffa.segments_overlap_bp([(1, 10), (20, 30)], [(5, 25)]) == 6 + 6
    assert ffa.segments_overlap_bp([(1, 10)], [(11, 20)]) == 0


def test_valid_orf_tristate(ffa):
    assert ffa._valid_orf("valid_ORF=True") is True
    assert ffa._valid_orf("valid_ORFs=0") is False
    assert ffa._valid_orf("ID=g1") is None


def test_find_cds_overlaps_uses_cds_not_gene_span(ffa):
    genes = {
        "a": _gene("a", 1, 500, [(1, 100), (400, 500)]),
        "b": _gene("b", 150, 350, [(150, 350)]),  # inside a's intron: gene spans overlap, CDS don't
        "c": _gene("c", 450, 700, [(450, 700)]),
        "d": _gene("d", 450, 700, [(450, 700)], seqname="chr2"),
    }
    assert ffa.find_cds_overlaps(genes, min_overlap_bp=20) == [("a", "c", 51)]
    assert ffa.find_cds_overlaps(genes, min_overlap_bp=60) == []


def test_classify_minor_position(ffa):
    bigger = _gene("big", 1, 1000, [(1, 100), (400, 500), (900, 1000)])
    assert ffa.classify_minor_position(bigger, [(950, 1100)]) == "flank"
    assert ffa.classify_minor_position(bigger, [(450, 460)]) == "inside"
    assert ffa.classify_minor_position(_gene("one", 1, 100, [(1, 100)]), [(90, 120)]) == "single_exon"


def test_is_qualifying_overlap_minor_flank_is_not_a_conflict(ffa):
    bigger = _gene("big", 1, 1000, [(1, 100), (400, 500), (900, 1000)])
    flank = _gene("s", 990, 1150, [(990, 1150)])  # 11 bp of a 161 bp CDS, on the big gene's last segment
    qualifies, pct = ffa.is_qualifying_overlap(bigger, flank, 11)
    assert qualifies is False
    assert pct == pytest.approx(11 / 161)

    contained = _gene("s2", 410, 480, [(410, 480)])  # fully inside an internal segment
    assert ffa.is_qualifying_overlap(bigger, contained, 71) == (True, 1.0)


def test_panel_check(ffa):
    present = {"rare": {"A", "B"}, "common": {"A", "B", "C", "D"}}
    assert ffa.panel_check("rare", "x", 0.5, present, {"rare": {"A"}}, 2, "") is not None
    assert ffa.panel_check("common", "x", 0.5, present, {"common": {"A"}}, 2, "") is None
    assert ffa.panel_check("common", "x", 0.5, present, {"common": {"A", "B"}}, 2, "") is not None


def test_assign_new_ids_orders_by_position(ffa):
    genes = {"g2": _gene("g2", 500, 600, []), "g1": _gene("g1", 100, 200, []), "g3": _gene("g3", 1, 50, [], "chr2")}
    assert ffa.assign_new_ids(genes, {"g1", "g2", "g3"}, 3, "Tg") == {
        "g1": "Tg3g00001",
        "g2": "Tg3g00002",
        "g3": "Tg3g00003",
    }


def test_parse_gff_picks_longest_mrna_and_reads_evidence(ffa, write_file):
    text = (
        gff_line("chr1", "protein_coding_gene", 1, 500, "+", "ID=g1;description=heat+shock%20protein;coverage=0.8")
        + gff_line("chr1", "mRNA", 1, 500, "+", "ID=m1;Parent=g1")
        + gff_line("chr1", "CDS", 1, 100, "+", "ID=c1;Parent=m1")
        + gff_line("chr1", "mRNA", 1, 500, "+", "ID=m2;Parent=g1")
        + gff_line("chr1", "CDS", 1, 100, "+", "ID=c2;Parent=m2")
        + gff_line("chr1", "CDS", 300, 500, "+", "ID=c3;Parent=m2;Dbxref=InterPro:IPR1")
        + gff_line("chr1", "protein_coding_gene", 600, 700, "+", "ID=g2;description=hypothetical protein")
    )
    genes = ffa.parse_gff(write_file("final.gff", text))
    g1 = genes["g1"]
    assert g1["cds_segments"] == [(1, 100), (300, 500)]
    assert g1["cds_bp"] == 301
    assert g1["has_domain"] is True
    assert g1["description"] == "heat shock protein"
    assert g1["coverage"] == pytest.approx(0.8)
    assert ffa.has_evidence(g1, require_domain_only=False)
    assert not ffa.has_evidence(genes["g2"], require_domain_only=False)  # generic description
    assert genes["g2"]["cds_bp"] == 0
