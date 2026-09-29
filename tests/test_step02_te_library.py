"""workflow/step02_te_annotation_clean_up/te_library - the TE-library pipeline's config and job wiring.

Run from a project folder, everything it writes goes under output/step02_te_annotation_clean_up/te_library/.
"""

import os
import subprocess

from conftest import REPO_ROOT

TE = "workflow/step02_te_annotation_clean_up/te_library"


def _project(tmp_path):
    project = tmp_path / "project"
    (project / "genomes").mkdir(parents=True)
    for genome in ("ME49", "GT1"):
        (project / "genomes" / f"{genome}.fna").write_text(f">{genome}\nACGT\n")
    return project


def test_prepare_genome_table_writes_under_step02_out(tmp_path):
    """The bash-run helpers find config.sh next to themselves; PROJECT_DIR is the folder they run from."""
    project = _project(tmp_path)
    env = {**os.environ, "PIPELINE_CONFIG": str(tmp_path / "none.sh")}
    env.pop("SLURM_SUBMIT_DIR", None)
    env.pop("PIPELINE_DIR", None)
    result = subprocess.run(
        ["bash", str(REPO_ROOT / TE / "00_prepare_genome_table.sh")],
        cwd=project,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    table = project / "output" / "step02_te_annotation_clean_up" / "te_library" / "config" / "genomes.tsv"
    assert sorted(table.read_text().splitlines()) == [
        f"GT1\t{project}/genomes/GT1.fna",
        f"ME49\t{project}/genomes/ME49.fna",
    ]
    assert not (REPO_ROOT / TE / "config").exists()  # nothing written into the checkout


def test_step_job_finds_config_through_pipeline_dir(run_job, tmp_path):
    """A spooled StepN job loads te_library/config.sh from the checkout and reads the project's outputs."""
    project = _project(tmp_path)
    result = run_job(f"{TE}/Step1_run_repeatmodeler.slurm.sh", project)
    assert result.returncode == 1
    table = project / "output" / "step02_te_annotation_clean_up" / "te_library" / "config" / "genomes.tsv"
    assert f"genome table not found: {table}" in result.stderr
