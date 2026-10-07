"""lib/pangenome_common - helpers shared across step scripts."""

from pathlib import Path

import pytest

from pangenome_common.gff import (
    biotype_category,
    find_gff,
    gene_description,
    gff_seqids,
    infer_biotype,
    parse_attributes,
)
from pangenome_common.intervals import merged_len, segments_overlap_bp
from pangenome_common.io import read_id_list, read_paf, read_strains


def test_merged_len_half_open():
    assert merged_len([]) == 0
    assert merged_len([(0, 10)]) == 10
    assert merged_len([(20, 30), (0, 10), (5, 15)]) == 25  # [0,15) + [20,30)
    assert merged_len([(0, 10), (10, 20)]) == 20  # touching half-open intervals don't double count
    assert merged_len([(0, 100), (10, 20)]) == 100  # nested


def test_segments_overlap_bp_closed():
    assert segments_overlap_bp([(1, 10)], [(10, 20)]) == 1  # 1-based inclusive: position 10 is shared
    assert segments_overlap_bp([(1, 10), (20, 30)], [(5, 25)]) == 12
    assert segments_overlap_bp([(1, 10)], [(11, 20)]) == 0
    assert segments_overlap_bp([], [(1, 5)]) == 0


def test_read_strains_from_file_and_string(tmp_path):
    order = tmp_path / "genome_order.txt"
    order.write_text("# reference first\nREF\n\n  S2  \nS3\n")
    assert read_strains(str(order)) == ["REF", "S2", "S3"]
    assert read_strains("REF, S2,,S3") == ["REF", "S2", "S3"]


def test_find_gff_extension_order_and_return_type(tmp_path):
    (tmp_path / "G1.gff").write_text("")
    (tmp_path / "G1.gff3").write_text("")
    assert find_gff(str(tmp_path), "G1") == str(tmp_path / "G1.gff3")  # default prefers .gff3
    assert find_gff(tmp_path, "G1", exts=(".gff", ".gff3")) == tmp_path / "G1.gff"
    assert isinstance(find_gff(tmp_path, "G1"), Path)
    assert find_gff(str(tmp_path), "missing") is None


def test_find_gff_optional_gz_and_case_insensitive(tmp_path):
    (tmp_path / "Af293.gff3.gz").write_text("")
    (tmp_path / "FungiDB-68_AF293.gff").write_text("")
    assert find_gff(str(tmp_path), "Af293") is None
    assert find_gff(str(tmp_path), "Af293", exts=(".gff3.gz",)) == str(tmp_path / "Af293.gff3.gz")
    assert find_gff(tmp_path, "fungidb-68_af293", case_insensitive=True) == tmp_path / "FungiDB-68_AF293.gff"


def test_gff_seqids_column_one_only(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##sequence-region chrX 1 5\nchr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=g\n\nchr2\tsrc\tgene\t1\t5\t.\t+\t.\tID=h\n"
    )
    assert gff_seqids(gff) == {"chr1", "chr2"}


def test_read_id_list_takes_first_column(tmp_path):
    ids = tmp_path / "ids.txt"
    ids.write_text("# comment\ngeneA  some note\n\ngeneB\n")
    assert read_id_list(str(ids)) == ["geneA", "geneB"]
    assert read_id_list("geneA,geneB") == ["geneA", "geneB"]


def test_read_paf_filters_short_alignments(tmp_path):
    paf = tmp_path / "a.paf"
    paf.write_text(
        "ctg1\t1000\t0\t900\t+\tchr1\t5000\t100\t1000\t880\t900\t60\n"
        "ctg2\t1000\t0\t50\t-\tchr1\t5000\t100\t150\t50\t50\t60\n"
        "short line\n"
    )
    rows = read_paf(paf, min_aln=100)
    assert rows == [
        {
            "q": "ctg1",
            "qlen": 1000,
            "qs": 0,
            "qe": 900,
            "strand": "+",
            "t": "chr1",
            "tlen": 5000,
            "ts": 100,
            "te": 1000,
            "nmatch": 880,
            "alen": 900,
        }
    ]


def test_infer_biotype_attribute_wins_over_feature_type():
    assert infer_biotype("gene", {"gene_biotype": "tRNA"}) == "tRNA"
    assert infer_biotype("gene", {}) == "protein_coding"
    assert infer_biotype("ncRNA_gene", {}) == "ncRNA"
    assert infer_biotype("mystery", {}) == "NA"


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("protein_coding", "protein_coding"),
        ("mRNA", "protein_coding"),
        ('"pseudogene"', "pseudogene"),
        ("processed_pseudogene", "pseudogene"),
        ("rRNA_pseudogene", "pseudogene"),  # pseudogene checked before rRNA
        ("Mt_tRNA", "tRNA"),
        ("lncRNA", "ncRNA"),
        ("scaRNA", "ncRNA"),
        ("IG_C_gene", None),
    ],
)
def test_biotype_category(value, expected):
    assert biotype_category(value) == expected


def test_biotype_category_other_rna_label():
    assert biotype_category("snoRNA", other_rna_label="ncRNA_other") == "ncRNA_other"
    assert biotype_category("tRNA", other_rna_label="ncRNA_other") == "tRNA"


def test_parse_attributes():
    assert parse_attributes("ID=g1;Note=a=b;flag;Name=x \r") == {"ID": "g1", "Note": "a=b", "Name": "x"}


def test_gene_description_first_non_empty_and_decoded():
    assert gene_description({"description": "hypothetical+protein"}) == "hypothetical protein"
    assert gene_description({"description": "heat%20shock%2Bcofactor"}) == "heat shock+cofactor"
    assert gene_description({"description": "", "product": "kinase"}) == "kinase"
    assert gene_description({"Name": "abc"}) == "abc"
    assert gene_description({}) == ""
