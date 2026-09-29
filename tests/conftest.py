"""Shared test helpers.

load_script() imports pipeline code by its repo path: modules in src/pangenome_pipeline by their
package name, and the loose scripts that remain outside it (species_specific/, scripts/) by
file path, with their folder on sys.path for sibling imports. run_job() runs a workflow/ job script
the way sbatch does.
"""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

REPO_ROOT = Path(__file__).resolve().parent.parent


def module_name(relpath):
    """src/pangenome_pipeline/x/y.py -> pangenome_pipeline.x.y (None for files outside src/)."""
    if not relpath.startswith("src/") or not relpath.endswith(".py"):
        return None
    return relpath[len("src/") : -len(".py")].replace("/", ".")


def _load(relpath):
    if module_name(relpath):
        return importlib.import_module(module_name(relpath))
    path = REPO_ROOT / relpath
    script_dir = str(path.parent)
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)  # sibling imports, e.g. `from plot_pancore import read_presence`
    name = "pipeline_" + path.stem
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def load_script():
    return _load


@pytest.fixture
def write_file(tmp_path):
    """Write text to tmp_path/<name> and return the path."""

    def _write(name, text):
        path = tmp_path / name
        path.write_text(text)
        return path

    return _write


def gff_line(seqid, ftype, start, end, strand, attrs, source="test"):
    return "\t".join([seqid, source, ftype, str(start), str(end), ".", strand, ".", attrs]) + "\n"


@pytest.fixture
def run_job(tmp_path):
    """Run a workflow/ job script the way sbatch does: a spool copy, submitted from a project folder.

    run_job("workflow/stepNN_x/job.sh", project, *args, **env) -> CompletedProcess. The pipeline is
    found through PIPELINE_DIR; `module` and micromamba activation are no-ops, and python3 is this
    interpreter. The checkout's own pipeline_config.sh (if any) is not loaded.
    """
    stubs = tmp_path / "stub_bin"
    stubs.mkdir()
    (stubs / "module").write_text("#!/bin/sh\nexit 0\n")
    (stubs / "micromamba").write_text("#!/bin/sh\necho 'micromamba() { :; }'\n")
    for stub in ("module", "micromamba"):
        (stubs / stub).chmod(0o755)
    (stubs / "python3").symlink_to(sys.executable)

    def _run(script, project, *args, **env):
        spool = tmp_path / "spool"
        spool.mkdir(exist_ok=True)
        job = spool / "slurm_script"
        job.write_text((REPO_ROOT / script).read_text())
        full_env = {
            **os.environ,
            "PATH": f"{stubs}:{os.environ['PATH']}",
            "PIPELINE_DIR": str(REPO_ROOT),
            "SLURM_SUBMIT_DIR": str(project),
            "PIPELINE_CONFIG": str(tmp_path / "no_pipeline_config.sh"),
            "MAMBA_EXE": str(stubs / "micromamba"),
            "MAMBA_ROOT_PREFIX": str(tmp_path),
            **env,
        }
        return subprocess.run(
            ["bash", str(job), *args], cwd=project, env=full_env, capture_output=True, text=True, check=False
        )

    return _run
