"""lib/slurm_common.sh - the helpers every job script sources (find_script, activate_env, config)."""

import os
import subprocess

import pytest
from conftest import REPO_ROOT

COMMON = REPO_ROOT / "lib" / "slurm_common.sh"

# How every job script finds the pipeline (see lib/slurm_common.sh). tests/test_workflow_layout.py
# checks that each job under workflow/ starts with it.
PREAMBLE = (
    '_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"\n'
    '[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; '
    'while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }\n'
    'source "$_pd/lib/slurm_common.sh" || '
    '{ echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }\n'
)


def _bash(script, cwd, **env):
    return subprocess.run(
        ["bash", "-c", script], cwd=cwd, env={**os.environ, **env}, capture_output=True, text=True, check=False
    )


@pytest.fixture
def project(tmp_path):
    """A species project folder: inputs only, no pipeline code in it."""
    project = tmp_path / "project"
    project.mkdir()
    (project / "genome_order.txt").write_text("ME49\n")
    return project


def _job(folder, body):
    """A job script at folder/slurm_script that starts with the standard preamble."""
    folder.mkdir(parents=True, exist_ok=True)
    job = folder / "slurm_script"
    job.write_text(f"#!/bin/bash\nset -euo pipefail\n{PREAMBLE}{body}\n")
    return job


def _run_job(job, cwd, **env):
    env = {k: v for k, v in {**os.environ, **env}.items() if v is not None}
    return subprocess.run(["bash", str(job)], cwd=cwd, env=env, capture_output=True, text=True, check=False)


def test_job_finds_pipeline_from_pipeline_dir(project, tmp_path):
    """Under sbatch the job runs from a spool copy, so it needs PIPELINE_DIR (passed through by sbatch)."""
    job = _job(tmp_path / "spool", 'echo "$PIPELINE_DIR|$PROJECT_DIR|$OUTPUT_DIR|$STEP05_OUT"')
    result = _run_job(job, project, PIPELINE_DIR=str(REPO_ROOT), SLURM_SUBMIT_DIR=str(project))
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == (
        f"{REPO_ROOT}|{project}|{project}/output|{project}/output/step05_cross_genome_annotation_integration"
    )


def test_job_finds_pipeline_from_its_own_location(project, tmp_path):
    """Run directly with bash, a job under workflow/stepNN_*/ finds the checkout it sits in."""
    checkout = tmp_path / "checkout"
    (checkout / "lib").mkdir(parents=True)
    (checkout / "lib" / "slurm_common.sh").write_text(COMMON.read_text())
    job = _job(checkout / "workflow" / "step01_x", 'echo "$PIPELINE_DIR"')
    result = _run_job(job, project, PIPELINE_DIR=None, SLURM_SUBMIT_DIR=None)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(checkout)


def test_job_finds_pipeline_above_the_submit_folder(tmp_path):
    """A checkout cloned into the project root is still found from the submit folder."""
    project = tmp_path / "clone"
    (project / "lib").mkdir(parents=True)
    (project / "lib" / "slurm_common.sh").symlink_to(REPO_ROOT / "lib" / "slurm_common.sh")
    job = _job(tmp_path / "spool", 'echo "$PROJECT_DIR"')
    result = _run_job(job, project, PIPELINE_DIR=None, SLURM_SUBMIT_DIR=str(project))
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(project)


def test_job_without_pipeline_stops(tmp_path):
    job = _job(tmp_path / "spool", "echo still-running")
    result = _run_job(job, tmp_path, PIPELINE_DIR=None, SLURM_SUBMIT_DIR=str(tmp_path))
    assert result.returncode == 1
    assert "export PIPELINE_DIR" in result.stderr
    assert "still-running" not in result.stdout


def test_output_dir_setting(project):
    script = f'source {COMMON}; echo "$OUTPUT_DIR|$STEP01_OUT|$STEP08_OUT|$ALT_PANLOCI_OUT"'
    relative = _bash(script, project, OUTPUT_DIR="run1").stdout.strip()
    assert relative == (
        f"{project}/run1|{project}/run1/step01_genome_selection_qc|"
        f"{project}/run1/step08_pangenome_visualisation|{project}/run1/alternative_get_panloci"
    )
    absolute = _bash(script, project, OUTPUT_DIR="/scratch/x", STEP08_OUT="/elsewhere").stdout.strip()
    assert absolute == "/scratch/x|/scratch/x/step01_genome_selection_qc|/elsewhere|/scratch/x/alternative_get_panloci"


def test_project_config_overrides_pipeline_config(project, tmp_path):
    pipeline_cfg = tmp_path / "cfg.sh"
    pipeline_cfg.write_text("GET_PANGENES_ENV=/from/pipeline\nENV_SYNTENY=pipeline_synteny\n")
    (project / "pipeline_config.sh").write_text("ENV_SYNTENY=project_synteny\nOUTPUT_DIR=/proj/out\n")
    script = f'source {COMMON}; echo "$GET_PANGENES_ENV $ENV_SYNTENY $STEP04_OUT"'
    result = _bash(script, project, PIPELINE_CONFIG=str(pipeline_cfg))
    assert result.stdout.strip() == "/from/pipeline project_synteny /proj/out/step04_interproscan"


def test_find_script_resolves_from_the_pipeline_root(project):
    result = _bash(f'source {COMMON}; find_script F lib/slurm_common.sh; echo "$F"', project)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f"{REPO_ROOT}/lib/slurm_common.sh"


def test_find_script_missing_stops_the_job(project):
    result = _bash(f"set -u; source {COMMON}; find_script PY workflow/nope.py; echo still-running", project)
    assert result.returncode == 1
    assert "cannot find workflow/nope.py in the pipeline" in result.stderr
    assert "still-running" not in result.stdout


def test_require_config_names_missing_settings(project):
    result = _bash(f"source {COMMON}; FOO=1; require_config FOO BAR BAZ", project, BAR="", BAZ="")
    assert result.returncode == 1
    assert "BAR BAZ not set" in result.stderr
    assert "pipeline_config.example.sh" in result.stderr


def test_config_file_is_loaded_and_env_can_override(project, tmp_path):
    config = tmp_path / "cfg.sh"
    config.write_text('GET_PANGENES_ENV="${GET_PANGENES_ENV:-/from/config}"\nENV_SYNTENY=my_synteny\n')
    script = f'source {COMMON}; echo "$GET_PANGENES_ENV $ENV_SYNTENY $ENV_BIOINFO"'
    assert _bash(script, project, PIPELINE_CONFIG=str(config)).stdout.strip() == "/from/config my_synteny bioinfo"
    overridden = _bash(script, project, PIPELINE_CONFIG=str(config), GET_PANGENES_ENV="/from/env").stdout.strip()
    assert overridden == "/from/env my_synteny bioinfo"


@pytest.fixture
def fake_micromamba(tmp_path):
    """A stand-in for micromamba whose shell hook defines a `micromamba activate` that records the env."""
    exe = tmp_path / "micromamba"
    exe.write_text(
        "#!/bin/bash\n"
        'echo \'micromamba() { [[ "$1" == activate ]] && [[ "$2" != missing ]] && export ACTIVE_ENV="$2"; }\'\n'
    )
    exe.chmod(0o755)
    return exe


@pytest.mark.parametrize("strict", ["set -euo pipefail", "set -eo pipefail"])
def test_activate_env_keeps_nounset_setting(project, fake_micromamba, tmp_path, strict):
    script = (
        f"{strict}; source {COMMON}; activate_env my_env; "
        'echo "$ACTIVE_ENV"; [[ $- == *u* ]] && echo nounset-on || echo nounset-off'
    )
    result = _bash(script, project, MAMBA_EXE=str(fake_micromamba), MAMBA_ROOT_PREFIX=str(tmp_path))
    assert result.returncode == 0, result.stderr
    expected = "nounset-on" if "u" in strict.split()[1] else "nounset-off"
    assert result.stdout.split() == ["my_env", expected]


def test_activate_env_errors(project, fake_micromamba, tmp_path):
    no_config = _bash(f"source {COMMON}; activate_env x", project, MAMBA_EXE="", MAMBA_ROOT_PREFIX="")
    assert no_config.returncode == 1 and "MAMBA_EXE MAMBA_ROOT_PREFIX not set" in no_config.stderr

    env = {"MAMBA_EXE": str(fake_micromamba), "MAMBA_ROOT_PREFIX": str(tmp_path)}
    missing = _bash(f"source {COMMON}; activate_env missing", project, **env)
    assert missing.returncode == 1 and "could not activate micromamba environment 'missing'" in missing.stderr
