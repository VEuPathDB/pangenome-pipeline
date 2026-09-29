"""step01_genome_selection_qc/drop_contigs.py - removing contigs from assemblies and GFFs."""

import io

import pytest
from conftest import gff_line

SCRIPT = "src/pangenome_pipeline/step01_genome_selection_qc/drop_contigs.py"


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


def test_drop_contigs_job_reads_the_list_from_step01_out(run_job, tmp_path):
    """Submitted from a project folder: the list is in output/step01_*/, genomes/ + gff_originals/ edited in place."""
    project = tmp_path / "project"
    (project / "genomes").mkdir(parents=True)
    (project / "gff_originals").mkdir()
    (project / "genomes" / "G1.fna").write_text(">chr1\nACGT\n>mito\nAAAA\n")
    (project / "gff_originals" / "G1.gff3").write_text(
        gff_line("chr1", "gene", 1, 4, "+", "ID=g1") + gff_line("mito", "gene", 1, 4, "+", "ID=g2")
    )
    step_out = project / "output" / "step01_genome_selection_qc"
    step_out.mkdir(parents=True)
    (step_out / "drop_contigs.txt").write_text("G1 mito\n")

    result = run_job("workflow/step01_genome_selection_qc/drop_contigs_slurm.sh", project)
    assert result.returncode == 0, result.stderr
    assert (project / "genomes" / "G1.fna").read_text() == ">chr1\nACGT\n"
    assert "g2" not in (project / "gff_originals" / "G1.gff3").read_text()


def test_drop_contigs_job_without_list_names_the_new_path(run_job, tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    result = run_job("workflow/step01_genome_selection_qc/drop_contigs_slurm.sh", project)
    assert result.returncode == 1
    assert f"{project}/output/step01_genome_selection_qc/drop_contigs.txt not found" in result.stderr
