"""workflow/step06_get_pangenes - the 7a/7b wrappers and their shared body (get_pangenes_common.sh)."""

import os
import subprocess

import pytest
from conftest import REPO_ROOT

JOBS = "workflow/step06_get_pangenes"
SCRIPTS = REPO_ROOT / JOBS
STEP05 = "output/step05_cross_genome_annotation_integration"
STEP06 = "output/step06_get_pangenes"


@pytest.fixture
def project(tmp_path):
    """A minimal project folder with two genomes' step-05 outputs (in output/step05_*/) and assemblies."""
    project = tmp_path / "project"
    (project / "genomes").mkdir(parents=True)
    (project / "genome_order.txt").write_text("# reference first\nME49\n\nGT1\n")
    for sub, suffix in (
        ("filtered_annotation", "_final_annotation_filtered.gff"),
        ("tagged_annotation", "_final_annotation_tagged.gff"),
    ):
        (project / STEP05 / sub).mkdir(parents=True)
        for genome in ("ME49", "GT1"):
            (project / STEP05 / sub / f"{genome}{suffix}").write_text(f"{sub} {genome}\n")
    for genome in ("ME49", "GT1"):
        (project / "genomes" / f"{genome}.fna").write_text(f">{genome}\nACGT\n")
    return project


def _bash(args, cwd, **env):
    return subprocess.run(
        ["bash", *args], cwd=cwd, env={**os.environ, **env}, capture_output=True, text=True, check=False
    )


@pytest.mark.parametrize(("variant", "annotation_dir"), [("7a", "filtered_annotation"), ("7b", "tagged_annotation")])
def test_link_job_symlinks_gff_and_fasta(run_job, project, variant, annotation_dir):
    result = run_job(f"{JOBS}/link_for_get_pangenes_{variant}.slurm.sh", project)
    assert result.returncode == 0, result.stderr
    assert "2 genome(s) linked" in result.stdout

    out = project / STEP06 / f"for_get_pangenes_{variant}"
    assert sorted(p.name for p in out.iterdir()) == [
        "GT1_for_get_pangenes.fna",
        "GT1_for_get_pangenes.gff",
        "ME49_for_get_pangenes.fna",
        "ME49_for_get_pangenes.gff",
    ]
    gff = out / "GT1_for_get_pangenes.gff"
    assert gff.is_symlink() and gff.read_text() == f"{annotation_dir} GT1\n"


def test_link_job_fails_on_missing_input(run_job, project):
    (project / "genomes" / "GT1.fna").unlink()
    result = run_job(f"{JOBS}/link_for_get_pangenes_7a.slurm.sh", project)
    assert result.returncode != 0
    assert "missing genomes/GT1.fna" in result.stderr


@pytest.mark.parametrize(
    ("taxon_regex", "expected_args"),
    [
        ("", "-d in -t 0 -n 4 -B {env}/bedtools"),
        ("GT", "-d in -t 0 -s GT -n 4 -B {env}/bedtools"),
    ],
)
def test_run_get_pangenes_calls_tool_with_settings(tmp_path, taxon_regex, expected_args):
    env_bin = tmp_path / "env_bin"
    env_bin.mkdir()
    fake = env_bin / "get_pangenes.pl"
    fake.write_text('#!/bin/bash\necho "$@" > "$(dirname "$0")/called_with"; pwd > "$(dirname "$0")/ran_in"\n')
    fake.chmod(0o755)
    (tmp_path / "out" / "in").mkdir(parents=True)
    script = (
        "set -euo pipefail\n"
        "module() { :; }\n"  # no environment modules outside the cluster
        f"source {REPO_ROOT / 'lib' / 'slurm_common.sh'}\n"
        f"source {SCRIPTS / 'get_pangenes_common.sh'}\n"
        f"GET_PANGENES_ENV={env_bin} INPUT_DIR={tmp_path}/out/in N_CPUS=4 TAXA_THRESHOLD=0 VARIANT=7a\n"
        "run_get_pangenes\n"
    )
    result = _bash(["-c", script], cwd=tmp_path, TAXON_REGEX=taxon_regex)
    assert result.returncode == 0, result.stderr
    assert "get_pangenes (7a) completed" in result.stdout
    assert (env_bin / "called_with").read_text().strip() == expected_args.format(env=env_bin)
    # run from INPUT_DIR's parent, so get_pangenes' <input>_pangenes/ lands next to INPUT_DIR
    assert (env_bin / "ran_in").read_text().strip() == str(tmp_path / "out")


def test_run_get_pangenes_without_env_uses_activate_env(tmp_path):
    """No GET_PANGENES_ENV: get_pangenes.pl comes from activate_env (a container, or PATH), without -B."""
    tools = tmp_path / "tools"
    tools.mkdir()
    fake = tools / "get_pangenes.pl"
    fake.write_text('#!/bin/bash\necho "$@" > "$(dirname "$0")/called_with"\n')
    fake.chmod(0o755)
    (tmp_path / "out" / "in").mkdir(parents=True)
    script = (
        "set -euo pipefail\n"
        "module() { :; }\n"
        f"source {REPO_ROOT / 'lib' / 'slurm_common.sh'}\n"
        f"source {SCRIPTS / 'get_pangenes_common.sh'}\n"
        f"INPUT_DIR={tmp_path}/out/in N_CPUS=4 TAXA_THRESHOLD=0 VARIANT=7a\n"
        "run_get_pangenes\n"
    )
    env = {"GET_PANGENES_ENV": "", "MAMBA_EXE": "", "CONDA_EXE": "", "CONTAINER_DIR": "", "ENV_BACKEND": "path"}
    result = _bash(["-c", script], cwd=tmp_path, PATH=f"{tools}:{os.environ['PATH']}", **env)
    assert result.returncode == 0, result.stderr
    assert "activate_env get_pangenes: path" in result.stderr
    assert (tools / "called_with").read_text().strip() == "-d in -t 0 -n 4"


def test_run_get_pangenes_says_how_to_provide_it(tmp_path):
    script = (
        "module() { :; }\n"
        f"source {REPO_ROOT / 'lib' / 'slurm_common.sh'}\n"
        f"source {SCRIPTS / 'get_pangenes_common.sh'}\n"
        "INPUT_DIR=in N_CPUS=4 TAXA_THRESHOLD=0 VARIANT=7a\n"
        "run_get_pangenes\n"
    )
    result = _bash(
        ["-c", script], cwd=tmp_path, GET_PANGENES_ENV="", ENV_BACKEND="path", PIPELINE_CONFIG=str(tmp_path / "none.sh")
    )
    assert result.returncode == 1
    assert "get_pangenes.pl not found: set GET_PANGENES_ENV" in result.stderr
