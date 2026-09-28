"""lib/slurm_common.sh - the helpers every job script sources (find_script, activate_env, config)."""

import os
import subprocess

import pytest
from conftest import REPO_ROOT

COMMON = REPO_ROOT / "lib" / "slurm_common.sh"


def _bash(script, cwd, **env):
    return subprocess.run(
        ["bash", "-c", script], cwd=cwd, env={**os.environ, **env}, capture_output=True, text=True, check=False
    )


@pytest.fixture
def project(tmp_path):
    """A project folder laid out like a real one: the pipeline's lib/ and a step folder in it."""
    (tmp_path / "lib").symlink_to(REPO_ROOT / "lib")
    (tmp_path / "src").symlink_to(REPO_ROOT / "src")
    step = tmp_path / "05_step" / "scripts"
    step.mkdir(parents=True)
    (step / "tool.py").write_text("")
    return tmp_path


def test_find_script_from_submit_dir(project, tmp_path_factory):
    """Under sbatch the wrapper runs from a spool copy, so the path comes from SLURM_SUBMIT_DIR."""
    spool = tmp_path_factory.mktemp("spool")
    job = spool / "slurm_script"
    job.write_text(
        '#!/bin/bash\nset -euo pipefail\nsource "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1\n'
        'find_script PY 05_step/scripts/tool.py\necho "$PY"\n'
    )
    result = subprocess.run(
        ["bash", str(job)],
        cwd=spool,
        env={**os.environ, "SLURM_SUBMIT_DIR": str(project)},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f"{project}/05_step/scripts/tool.py"


def test_find_script_missing_stops_the_job(project):
    result = _bash(f"set -u; source {COMMON}; find_script PY 05_step/scripts/nope.py; echo still-running", project)
    assert result.returncode == 1
    assert "cannot find nope.py (submit from the project folder)" in result.stderr
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
