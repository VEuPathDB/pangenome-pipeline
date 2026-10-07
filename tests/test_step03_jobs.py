"""03_denovo_Genome_Annotation jobs, submitted the way its README says: from inside the step folder.

They find the pipeline's lib/slurm_common.sh by looking upwards from the submit folder.
"""

import gzip
import os
import subprocess
import sys

import pytest
from conftest import REPO_ROOT

STEP = "03_denovo_Genome_Annotation"


@pytest.fixture
def step_dir(tmp_path):
    """<project>/03_denovo_Genome_Annotation with its scripts/, and the pipeline's lib/ at the project root."""
    project = tmp_path / "project"
    (project / STEP).mkdir(parents=True)
    (project / "lib").symlink_to(REPO_ROOT / "lib")
    (project / STEP / "scripts").symlink_to(REPO_ROOT / STEP / "scripts")
    stubs = tmp_path / "bin"
    stubs.mkdir()
    (stubs / "module").write_text("#!/bin/sh\nexit 0\n")
    (stubs / "module").chmod(0o755)
    (stubs / "python3").symlink_to(sys.executable)
    return project / STEP, stubs


def _submit(script, cwd, stubs, **env):
    return subprocess.run(
        ["bash", f"scripts/{script}"],
        cwd=cwd,
        env={**os.environ, "PATH": f"{stubs}:{os.environ['PATH']}", "SLURM_SUBMIT_DIR": str(cwd), **env},
        capture_output=True,
        text=True,
        check=False,
    )


def test_fix_headers_job_from_the_step_folder(step_dir):
    cwd, stubs = step_dir
    for side in ("left", "right"):
        with gzip.open(cwd / f"combined_{side}_trinity.fastq.gz", "wt") as fh:
            fh.write("@read1 extra words\nACGT\n+read1 extra words\nIIII\n")
    result = _submit("fix_headers.slurm.sh", cwd, stubs)
    assert result.returncode == 0, result.stderr
    with gzip.open(cwd / "combined_left_trinity.fixed_headers.fastq.gz", "rt") as fh:
        assert fh.read() == "@read1/1\nACGT\n+\nIIII\n"


def test_funannotate_job_finds_lib_from_the_step_folder(step_dir, tmp_path):
    """It gets past loading lib/slurm_common.sh and stops at the missing setting, not at the source line."""
    cwd, stubs = step_dir
    result = _submit(
        "busco.slurm.sh", cwd, stubs, FUNANNOTATE_SHARED="", PIPELINE_CONFIG=str(tmp_path / "no_config.sh")
    )
    assert result.returncode == 1
    assert "FUNANNOTATE_SHARED not set" in result.stderr
    assert "slurm_common.sh not found" not in result.stderr
