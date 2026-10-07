"""Baseline smoke tests: every pipeline script at least compiles, parses its CLI and has valid shell syntax.

These don't check what the scripts compute - that needs small GFF/FASTA fixtures - but they catch
broken imports (including the sibling-script imports), syntax errors and argparse mistakes.
"""

import os
import py_compile
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

# results/ only holds snapshot copies of pipeline scripts.
EXCLUDED_PREFIXES = ("results/",)


def tracked_files(*patterns):
    out = subprocess.run(
        ["git", "ls-files", *patterns], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    return sorted(f for f in out if not f.startswith(EXCLUDED_PREFIXES) and f not in ("noxfile.py",))


PY_FILES = [f for f in tracked_files("*.py") if not f.startswith("tests/")]
SHELL_FILES = tracked_files("*.sh", "*.slurm")
ARGPARSE_FILES = [f for f in PY_FILES if "argparse" in (REPO_ROOT / f).read_text()]


@pytest.mark.parametrize("path", PY_FILES)
def test_python_compiles(path, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        py_compile.compile(str(REPO_ROOT / path), cfile=str(tmp_path / "out.pyc"), doraise=True)


@pytest.mark.parametrize("path", ARGPARSE_FILES)
def test_python_help(path):
    script = REPO_ROOT / path
    env = {**os.environ, "MPLBACKEND": "Agg", "PYTHONDONTWRITEBYTECODE": "1"}
    result = subprocess.run(
        [sys.executable, script.name, "--help"],
        cwd=script.parent,
        env=env,
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("path", SHELL_FILES)
def test_shell_syntax(path):
    result = subprocess.run(["bash", "-n", str(REPO_ROOT / path)], capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr


def test_uv_lock_matches_pyproject():
    """Jobs run the Python code from a venv made with `uv sync --frozen` - uv.lock must be current."""
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip("uv not installed")
    result = subprocess.run(
        [uv, "lock", "--check", "--offline"], cwd=REPO_ROOT, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, "uv.lock is out of date - run `uv lock`:\n" + result.stderr
