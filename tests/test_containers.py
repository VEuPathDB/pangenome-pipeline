"""containers/ - every tool environment the jobs activate can be provided as a container."""

import re
import subprocess

from conftest import REPO_ROOT

PULL = REPO_ROOT / "containers" / "pull_containers.sh"


def _activated_envs():
    files = subprocess.run(
        ["git", "ls-files", "*.sh", "*.slurm"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    envs = set()
    for f in files:
        envs |= set(re.findall(r'activate_env "\$(ENV_[A-Z_]+)"', (REPO_ROOT / f).read_text(errors="replace")))
    assert "ENV_SYNTENY" in envs, envs  # the search found the jobs
    return envs


def test_every_activated_env_has_a_container():
    pulled = set(re.findall(r'"\$(ENV_[A-Z_]+)[/.]', PULL.read_text()))
    built = {"ENV_MCHELPER"}  # containers/mchelper.def
    missing = _activated_envs() - pulled - built
    assert not missing, f"add an image for {sorted(missing)} to containers/pull_containers.sh"


def test_pull_containers_needs_a_folder(tmp_path):
    result = subprocess.run(
        ["bash", str(PULL)],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "PIPELINE_CONFIG": str(tmp_path / "none.sh")},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1 and "give the folder" in result.stderr
