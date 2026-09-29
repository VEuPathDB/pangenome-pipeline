"""workflow/step03_denovo_genome_annotation jobs, submitted from the project folder.

They work in $STEP03_OUT (output/step03_denovo_genome_annotation/): the step's own inputs are put there
and every output lands there.
"""

import gzip

import pytest

STEP = "workflow/step03_denovo_genome_annotation"


@pytest.fixture
def project(tmp_path):
    """A project folder whose step-03 working folder holds the combined RNA-seq reads."""
    project = tmp_path / "project"
    work = project / "output" / "step03_denovo_genome_annotation"
    work.mkdir(parents=True)
    for side in ("left", "right"):
        with gzip.open(work / f"combined_{side}_trinity.fastq.gz", "wt") as fh:
            fh.write("@read1 extra words\nACGT\n+read1 extra words\nIIII\n")
    return project


def test_fix_headers_job_works_in_step03_out(run_job, project):
    result = run_job(f"{STEP}/fix_headers.slurm.sh", project)
    assert result.returncode == 0, result.stderr
    work = project / "output" / "step03_denovo_genome_annotation"
    with gzip.open(work / "combined_left_trinity.fixed_headers.fastq.gz", "rt") as fh:
        assert fh.read() == "@read1/1\nACGT\n+\nIIII\n"


def test_fix_headers_job_honours_output_dir(run_job, project, tmp_path):
    """OUTPUT_DIR moves the whole working folder."""
    elsewhere = tmp_path / "scratch"
    (project / "output").rename(elsewhere)
    result = run_job(f"{STEP}/fix_headers.slurm.sh", project, OUTPUT_DIR=str(elsewhere))
    assert result.returncode == 0, result.stderr
    assert (elsewhere / "step03_denovo_genome_annotation" / "combined_right_trinity.fixed_headers.fastq.gz").exists()


def test_funannotate_job_stops_at_missing_setting(run_job, project):
    """It gets past loading lib/slurm_common.sh and stops at the missing setting, not at the source line."""
    result = run_job(f"{STEP}/busco.slurm.sh", project, FUNANNOTATE_SHARED="")
    assert result.returncode == 1
    assert "FUNANNOTATE_SHARED not set" in result.stderr
    assert "pipeline not found" not in result.stderr
