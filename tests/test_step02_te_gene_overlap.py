"""02_TE_annotation_Clean_Up/scripts/te_gene_overlap.py - interval maths and TE/gene overlap classification."""

import pytest

SCRIPT = "02_TE_annotation_Clean_Up/scripts/te_gene_overlap.py"


@pytest.fixture(scope="module")
def te(load_script):
    return load_script(SCRIPT)


def test_merge_intervals_merges_overlapping_and_touching(te):
    assert te.merge_intervals([(20, 30), (1, 5), (6, 10), (25, 40), (50, 60)]) == [(1, 10), (20, 40), (50, 60)]
    assert te.merge_intervals([]) == []


def test_merged_track_overlap_bp_and_ranges(te):
    track = te.MergedTrack(te.merge_intervals([(10, 20), (40, 50)]))
    assert track.overlap_bp([(15, 45)]) == 6 + 6  # 15-20 and 40-45, 1-based inclusive
    assert track.overlap_bp([(21, 39)]) == 0
    assert track.overlap_ranges(15, 45) == [(15, 20), (40, 45)]


def test_parse_attributes_splits_on_first_equals(te):
    assert te.parse_attributes("ID=g1;Note=a=b;flag") == {"ID": "g1", "Note": "a=b"}


def test_infer_biotype_prefers_explicit_attribute(te):
    assert te.infer_biotype("gene", {"gene_biotype": "tRNA"}) == "tRNA"
    assert te.infer_biotype("gene", {}) == "protein_coding"
    assert te.infer_biotype("ncRNA_gene", {}) == "ncRNA"
    assert te.infer_biotype("weird_type", {}) == "NA"


def test_gene_description_order_and_url_decoding(te):
    assert te.gene_description({"product": "kinase", "description": "heat%20shock"}) == "heat shock"
    assert te.gene_description({"Name": "abc"}) == "abc"
    assert te.gene_description({}) == ""


@pytest.mark.parametrize(
    ("strand", "expected"),
    [("+", "5_prime"), ("-", "3_prime")],  # same genomic hit is the 3' end of a minus-strand gene
)
def test_classify_overlap_position_respects_strand(te, strand, expected):
    track = te.MergedTrack([(1, 100)])
    label, *_ = te.classify_overlap_position([(1, 1000)], strand, track, overlap_bp=100, transcript_bp=1000)
    assert label == expected


def test_classify_overlap_position_short_cases(te):
    track = te.MergedTrack([(1, 900)])
    assert te.classify_overlap_position([(1, 1000)], "+", track, 900, 1000)[0] == "spans_gene"
    assert te.classify_overlap_position([(1, 100)], "+", te.MergedTrack([(1, 10)]), 10, 100)[0] == "short_transcript"
