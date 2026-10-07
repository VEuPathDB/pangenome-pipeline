"""lib/slurm_common.sh - the helpers every job script sources (find_script, activate_env, config)."""

import os
import subprocess
from pathlib import Path

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
    micromamba = {"ENV_BACKEND": "micromamba", "PIPELINE_VENV": "/nonexistent"}
    no_config = _bash(f"source {COMMON}; activate_env x", project, MAMBA_EXE="", MAMBA_ROOT_PREFIX="", **micromamba)
    assert no_config.returncode == 1 and "MAMBA_EXE MAMBA_ROOT_PREFIX not set" in no_config.stderr

    env = {"MAMBA_EXE": str(fake_micromamba), "MAMBA_ROOT_PREFIX": str(tmp_path), **micromamba}
    missing = _bash(f"source {COMMON}; activate_env missing", project, **env)
    assert missing.returncode == 1 and "could not activate micromamba environment 'missing'" in missing.stderr

    bad = _bash(f"source {COMMON}; activate_env x", project, ENV_BACKEND="docker")
    assert bad.returncode == 1 and "ENV_BACKEND must be" in bad.stderr

    no_image = _bash(f"source {COMMON}; activate_env x", project, ENV_BACKEND="container", CONTAINER_DIR=str(tmp_path))
    assert no_image.returncode == 1 and "no container for 'x'" in no_image.stderr


def _no_envs(tmp_path):
    """Settings with no micromamba, conda or containers configured."""
    path = ":".join(d for d in os.environ["PATH"].split(":") if not (Path(d) / "conda").exists())
    return {"MAMBA_EXE": "", "CONDA_EXE": "", "CONTAINER_DIR": "", "PIPELINE_VENV": "/nonexistent", "PATH": path}


def test_activate_env_falls_back_to_path(project, tmp_path):
    """Nothing configured: the job uses what's already on PATH (modules etc.) and says so."""
    result = _bash(f"set -u; source {COMMON}; activate_env synteny", project, **_no_envs(tmp_path))
    assert result.returncode == 0, result.stderr
    assert "activate_env synteny: path" in result.stderr


def test_activate_env_puts_venv_python_on_path(project, tmp_path):
    venv = tmp_path / "venv"
    (venv / "bin").mkdir(parents=True)
    (venv / "bin" / "python3").write_text("#!/bin/bash\necho venv-python\n")
    (venv / "bin" / "python3").chmod(0o755)
    env = {**_no_envs(tmp_path), "PIPELINE_VENV": str(venv)}
    result = _bash(f"source {COMMON}; activate_env synteny; python3", project, **env)
    assert result.stdout.strip() == "venv-python", result.stderr


@pytest.fixture
def fake_apptainer(tmp_path):
    """A stand-in for apptainer: `exec IMG sh -c ...` lists the image's tools (the image file's lines);
    `exec --bind B IMG TOOL ARGS` prints what would run."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "apptainer"
    exe.write_text(
        "#!/bin/bash\n"
        '[[ "$1" == --silent ]] && shift; shift\n'
        'if [[ "$1" == --bind ]]; then echo "bind=$2 img=$(basename "$3") ${*:4}"; else cat "$1"; fi\n'
    )
    exe.chmod(0o755)
    return bindir


def test_activate_env_container(project, tmp_path, fake_apptainer):
    """A folder of images per environment: each image's tools get a wrapper, the first image wins."""
    images = tmp_path / "containers" / "synteny"
    images.mkdir(parents=True)
    (images / "a_minimap2.sif").write_text("minimap2\nsamtools\nbash\n")
    (images / "b_skani.sif").write_text("skani /opt/skani/skani\nsamtools\nbad|name\n")
    env = _no_envs(tmp_path)
    env.update(
        CONTAINER_DIR=str(tmp_path / "containers"),
        PATH=f"{fake_apptainer}:{env['PATH']}",
        TMPDIR=str(tmp_path),
        CONTAINER_BINDS="/data",
    )
    script = f"source {COMMON}; activate_env synteny; minimap2 -x asm10; samtools faidx; skani dist"
    result = _bash(script, project, **env)
    assert result.returncode == 0, result.stderr
    assert "activate_env synteny: container" in result.stderr
    shims = Path(
        _bash(f"source {COMMON}; activate_env synteny; dirname $(command -v minimap2)", project, **env).stdout.strip()
    )
    assert not (shims / "bash").exists(), "the host's own programs are not wrapped"
    binds = f"{project},{REPO_ROOT},/data"
    assert result.stdout.splitlines() == [
        f"bind={binds} img=a_minimap2.sif minimap2 -x asm10",
        f"bind={binds} img=a_minimap2.sif samtools faidx",
        f"bind={binds} img=b_skani.sif /opt/skani/skani dist",  # called at its real path
    ]
