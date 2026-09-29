"""species_specific/ jobs, submitted from a project folder like the pipeline's own."""

from conftest import gff_line

STEP01 = "output/step01_genome_selection_qc"


def test_tgondii_drop_short_contigs_writes_under_step01_out(run_job, tmp_path):
    project = tmp_path / "project"
    (project / "genomes").mkdir(parents=True)
    (project / "gff_originals").mkdir()
    (project / "genome_order.txt").write_text("G1\n")
    (project / "genomes" / "G1.fna").write_text(">chr1\nACGTACGTAC\n>tiny\nAC\n")
    (project / "genomes" / "G1.fna.fai").write_text("chr1\t10\t6\t10\t11\ntiny\t2\t23\t2\t3\n")  # no samtools needed
    (project / "gff_originals" / "G1.gff3").write_text(
        gff_line("chr1", "gene", 1, 10, "+", "ID=g1") + gff_line("tiny", "gene", 1, 2, "+", "ID=g2")
    )

    result = run_job("species_specific/Tgondii/step01_genome_selection_qc/drop_short_contigs.slurm.sh", project)
    assert result.returncode == 0, result.stderr
    step01 = project / STEP01
    # both contigs are under the script's MIN_LEN (10 kb)
    assert (step01 / "drop_short_contigs.txt").read_text().splitlines()[1:] == ["G1\tchr1", "G1\ttiny"]
    assert (step01 / "backup_before_short_contig_drop" / "genomes" / "G1.fna").read_text().startswith(">chr1")
    assert ">tiny" not in (project / "genomes" / "G1.fna").read_text()
    assert "summarize_dropped_genes.py not found - skipping the gene-cost summary" in result.stderr
