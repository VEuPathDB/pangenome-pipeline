"""step05_cross_genome_annotation_integration - Liftoff transfer filtering and overlap resolution."""

import os
import subprocess
import sys

import pytest
from conftest import gff_line

SCRIPTS = "src/pangenome_pipeline/step05_cross_genome_annotation_integration/"


@pytest.fixture(scope="module")
def good(load_script):
    return load_script(SCRIPTS + "get_good_transfer_ids.py")


@pytest.fixture(scope="module")
def ffa(load_script):
    return load_script(SCRIPTS + "filter_final_annotation.py")


@pytest.fixture(scope="module")
def common(load_script):
    return load_script(SCRIPTS + "annotation_common.py")


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
# annotation_common.py (shared by filter, tag and the QC plots)
#


def test_valid_orf_tristate(common):
    assert common.valid_orf("valid_ORF=True") is True
    assert common.valid_orf("valid_ORFs=0") is False
    assert common.valid_orf("ID=g1") is None


def test_attribute_readers_are_anchored(common):
    attrs = "sequence_ID=0.9;copy_num_ID=g1_0;ID=g1;Parent=m1;xextra_copy_number=7;extra_copy_number=2;xvalid_ORF=False;valid_ORF=True"
    assert common.ID_RE.search(attrs).group(1) == "g1"
    assert common.PARENT_RE.search(attrs).group(1) == "m1"
    assert common.int_attr(attrs, "extra_copy_number") == 2
    assert common.valid_orf(attrs) is True


@pytest.mark.parametrize(
    ("attrs", "expected"),
    [("coverage=0.95", 0.95), ("coverage=1e-05", 1e-05), ("coverage=1", 1.0), ("coverage=NA", None), ("ID=g", None)],
)
def test_float_attr_reads_whole_number(common, attrs, expected):
    assert common.float_attr(attrs, "coverage") == expected


@pytest.mark.parametrize(
    ("attrs", "key", "old", "new", "expected"),
    [
        ("ID=g1;Name=x", "ID", "g1", "Tst1g00001", "ID=Tst1g00001;Name=x"),
        ("copy_num_ID=g1;ID=g1", "ID", "g1", "N", "copy_num_ID=g1;ID=N"),  # not the copy_num_ID
        ("ID=g10;Note=ID=g1", "ID", "g1", "N", "ID=g10;Note=ID=g1"),  # no partial or nested match
        ("ID=m1;Parent=g1,g2", "Parent", "g1", "N", "ID=m1;Parent=N,g2"),
        ("ID=m1;Parent=g1", "Parent", "g1", "N", "ID=m1;Parent=N"),
    ],
)
def test_replace_attr_value(common, attrs, key, old, new, expected):
    assert common.replace_attr_value(attrs, key, old, new) == expected


def test_find_cds_overlaps_uses_cds_not_gene_span(common):
    genes = {
        "a": _gene("a", 1, 500, [(1, 100), (400, 500)]),
        "b": _gene("b", 150, 350, [(150, 350)]),  # inside a's intron: gene spans overlap, CDS don't
        "c": _gene("c", 450, 700, [(450, 700)]),
        "d": _gene("d", 450, 700, [(450, 700)], seqname="chr2"),
    }
    assert common.find_cds_overlaps(genes, min_overlap_bp=20) == [("a", "c", 51)]
    assert common.find_cds_overlaps(genes, min_overlap_bp=60) == []


def test_find_cds_overlaps_strand_flag_and_default_threshold(common):
    genes = {
        "a": _gene("a", 1, 100, [(1, 100)], strand="+"),
        "b": _gene("b", 100, 200, [(100, 200)], strand="-"),  # 1 bp overlap
    }
    assert common.find_cds_overlaps(genes) == [("a", "b", 1)]
    assert common.find_cds_overlaps(genes, with_strand=True) == [("a", "b", 1, False)]


def test_classify_minor_position(common):
    bigger = _gene("big", 1, 1000, [(1, 100), (400, 500), (900, 1000)])
    assert common.classify_minor_position(bigger, [(950, 1100)]) == "flank"
    assert common.classify_minor_position(bigger, [(450, 460)]) == "inside"
    assert common.classify_minor_position(_gene("one", 1, 100, [(1, 100)]), [(90, 120)]) == "single_exon"


#
# filter_final_annotation.py
#


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


def test_assign_new_ids_orders_by_position(common):
    genes = {"g2": _gene("g2", 500, 600, []), "g1": _gene("g1", 100, 200, []), "g3": _gene("g3", 1, 50, [], "chr2")}
    assert common.assign_new_ids(genes, {"g1", "g2", "g3"}, 3, "Tg") == {
        "g1": "Tg3g00001",
        "g2": "Tg3g00002",
        "g3": "Tg3g00003",
    }


def test_assign_new_ids_ties_are_reproducible(common):
    """Same contig and start: numbered by original ID, not by set/hash order."""
    script = (
        "import sys; sys.path.insert(0, sys.argv[1]); from annotation_common import assign_new_ids; "
        "genes = {g: {'seqname': 'chr1', 'start': 100} for g in ('geneA_copy1', 'geneB', 'geneC_copy2')}; "
        "print(sorted(assign_new_ids(genes, set(genes), 1, 'Afu').items()))"
    )
    outputs = {
        subprocess.run(
            [sys.executable, "-c", script, str(common.__file__).rsplit("/", 1)[0]],
            env={**os.environ, "PYTHONHASHSEED": str(seed)},
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        for seed in range(6)
    }
    assert outputs == {"[('geneA_copy1', 'Afu1g00001'), ('geneB', 'Afu1g00002'), ('geneC_copy2', 'Afu1g00003')]\n"}


def test_parse_gff_picks_longest_mrna_and_reads_evidence(common, ffa, write_file):
    text = (
        gff_line("chr1", "protein_coding_gene", 1, 500, "+", "ID=g1;description=heat+shock%20protein;coverage=0.8")
        + gff_line("chr1", "mRNA", 1, 500, "+", "ID=m1;Parent=g1")
        + gff_line("chr1", "CDS", 1, 100, "+", "ID=c1;Parent=m1")
        + gff_line("chr1", "mRNA", 1, 500, "+", "ID=m2;Parent=g1")
        + gff_line("chr1", "CDS", 1, 100, "+", "ID=c2;Parent=m2")
        + gff_line("chr1", "CDS", 300, 500, "+", "ID=c3;Parent=m2;Dbxref=InterPro:IPR1")
        + gff_line("chr1", "protein_coding_gene", 600, 700, "+", "ID=g2;description=hypothetical protein")
    )
    genes = common.parse_gff(write_file("final.gff", text))
    g1 = genes["g1"]
    assert g1["cds_segments"] == [(1, 100), (300, 500)]
    assert g1["cds_bp"] == 301
    assert g1["has_domain"] is True
    assert g1["description"] == "heat shock protein"
    assert g1["coverage"] == pytest.approx(0.8)
    assert ffa.has_evidence(g1, require_domain_only=False)
    assert not ffa.has_evidence(genes["g2"], require_domain_only=False)  # generic description
    assert genes["g2"]["cds_bp"] == 0
