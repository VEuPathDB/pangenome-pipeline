"""01_Genome_Selection_QC/scripts/drop_contigs.py - removing contigs from assemblies and GFFs."""

import io

import pytest
from conftest import gff_line

SCRIPT = "01_Genome_Selection_QC/scripts/drop_contigs.py"


@pytest.fixture(scope="module")
def dc(load_script):
    return load_script(SCRIPT)


def test_read_list_groups_contigs_by_genome_and_skips_comments(dc, write_file):
    path = write_file(
        "drop.txt",
        "# genome  contig\nAF293  mito_1\n\nGCA_1.1\tCP1.1  # trailing comment\nGCA_1.1 CP2.1\n",
    )
    assert dc.read_list(path) == {"AF293": {"mito_1"}, "GCA_1.1": {"CP1.1", "CP2.1"}}


def test_read_list_rejects_malformed_line(dc, write_file):
    path = write_file("drop.txt", "AF293 mito_1 extra\n")
    with pytest.raises(SystemExit):
        dc.read_list(path)


def test_filter_fasta_drops_named_records_and_reports_length(dc, write_file):
    path = write_file("g.fna", ">chr1 desc\nACGT\nACGT\n>mito\nAAAA\nAA\n>chr2\nGG\n")
    out = io.StringIO()
    dropped = dc.filter_fasta(path, {"mito"}, out)
    assert dropped == {"mito": 6}
    assert out.getvalue() == ">chr1 desc\nACGT\nACGT\n>chr2\nGG\n"


def test_filter_gff_drops_features_pragmas_and_embedded_fasta(dc, write_file):
    text = (
        "##gff-version 3\n"
        "##sequence-region chr1 1 100\n"
        "##sequence-region mito 1 50\n"
        + gff_line("chr1", "gene", 1, 10, "+", "ID=g1")
        + gff_line("mito", "gene", 1, 10, "+", "ID=g2")
        + gff_line("mito", "mRNA", 1, 10, "+", "ID=m2;Parent=g2")
        + "##FASTA\n>chr1\nACGT\n>mito\nTTTT\n"
    )
    path = write_file("g.gff3", text)
    out = io.StringIO()
    removed = dc.filter_gff(path, {"mito"}, out)
    result = out.getvalue()
    assert removed == 2
    assert "mito" not in result
    assert "ID=g1" in result and ">chr1\nACGT\n" in result


def test_gff_seqids_reads_columns_pragmas_and_fasta(dc, write_file):
    text = "##sequence-region chrA 1 5\n" + gff_line("chrB", "gene", 1, 5, "+", "ID=g") + "##FASTA\n>chrC\nAC\n"
    assert dc.gff_seqids(write_file("g.gff3", text)) == {"chrA", "chrB", "chrC"}
