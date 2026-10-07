"""setup_project.sh - creates the inputs and every log folder the jobs' #SBATCH lines write into."""

import re
import subprocess

from conftest import REPO_ROOT


def _job_log_dirs():
    files = subprocess.run(
        ["git", "ls-files", "*.sh", "*.slurm", ":!results"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    dirs = set()
    for f in files:
        text = (REPO_ROOT / f).read_text(errors="replace")
        dirs |= {d for d in re.findall(r"^#SBATCH --(?:output|error)=(\S+)/", text, re.MULTILINE)}
    assert "logs/slurm" in dirs, dirs  # the search found the jobs
    return dirs


def test_creates_inputs_and_every_log_folder(tmp_path):
    project = tmp_path / "project"
    result = subprocess.run(
        ["bash", str(REPO_ROOT / "setup_project.sh"), str(project)], capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
    for d in {"genomes", "gff_originals"} | _job_log_dirs():
        assert (project / d).is_dir(), d


def test_rerun_leaves_existing_files_alone(tmp_path):
    (tmp_path / "genomes").mkdir()
    (tmp_path / "genomes" / "ME49.fna").write_text(">c1\nACGT\n")
    for _ in range(2):
        result = subprocess.run(
            ["bash", str(REPO_ROOT / "setup_project.sh"), str(tmp_path)], capture_output=True, text=True, check=False
        )
        assert result.returncode == 0, result.stderr
    assert (tmp_path / "genomes" / "ME49.fna").read_text() == ">c1\nACGT\n"
