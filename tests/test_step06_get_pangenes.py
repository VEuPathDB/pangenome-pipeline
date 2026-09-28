"""06_GET_PANGENES - the 7a/7b wrappers and their shared body (get_pangenes_common.sh)."""

import os
import subprocess

import pytest
from conftest import REPO_ROOT

SCRIPTS = REPO_ROOT / "06_GET_PANGENES" / "scripts"


@pytest.fixture
def project(tmp_path):
    """A minimal project folder with two genomes' step-05 outputs and assemblies."""
    (tmp_path / "lib").symlink_to(REPO_ROOT / "lib")  # the pipeline is cloned into the project root
    (tmp_path / "genome_order.txt").write_text("# reference first\nME49\n\nGT1\n")
    for sub, suffix in (
        ("filtered_annotation", "_final_annotation_filtered.gff"),
        ("tagged_annotation", "_final_annotation_tagged.gff"),
    ):
        (tmp_path / sub).mkdir()
        for genome in ("ME49", "GT1"):
            (tmp_path / sub / f"{genome}{suffix}").write_text(f"{sub} {genome}\n")
    (tmp_path / "genomes").mkdir()
    for genome in ("ME49", "GT1"):
        (tmp_path / "genomes" / f"{genome}.fna").write_text(f">{genome}\nACGT\n")
    return tmp_path


def _bash(args, cwd, **env):
    return subprocess.run(
        ["bash", *args], cwd=cwd, env={**os.environ, **env}, capture_output=True, text=True, check=False
    )


@pytest.mark.parametrize(("variant", "annotation_dir"), [("7a", "filtered_annotation"), ("7b", "tagged_annotation")])
def test_link_wrapper_symlinks_gff_and_fasta(project, variant, annotation_dir):
    result = _bash([str(SCRIPTS / f"link_for_get_pangenes_{variant}.slurm.sh")], cwd=project)
    assert result.returncode == 0, result.stderr
    assert "2 genome(s) linked" in result.stdout

    out = project / f"for_get_pangenes_{variant}"
    assert sorted(p.name for p in out.iterdir()) == [
        "GT1_for_get_pangenes.fna",
        "GT1_for_get_pangenes.gff",
        "ME49_for_get_pangenes.fna",
        "ME49_for_get_pangenes.gff",
    ]
    gff = out / "GT1_for_get_pangenes.gff"
    assert gff.is_symlink() and gff.read_text() == f"{annotation_dir} GT1\n"


def test_link_wrapper_fails_on_missing_input(project):
    (project / "genomes" / "GT1.fna").unlink()
    result = _bash([str(SCRIPTS / "link_for_get_pangenes_7a.slurm.sh")], cwd=project)
    assert result.returncode != 0
    assert "missing genomes/GT1.fna" in result.stderr


def test_link_wrapper_finds_common_file_from_submit_dir(project, tmp_path_factory):
    """Under sbatch the wrapper runs from a spool copy - the shared file is found via SLURM_SUBMIT_DIR."""
    spool = tmp_path_factory.mktemp("spool")
    job = spool / "slurm_script"
    job.write_text((SCRIPTS / "link_for_get_pangenes_7a.slurm.sh").read_text())
    (project / "06_GET_PANGENES").symlink_to(REPO_ROOT / "06_GET_PANGENES")
    result = _bash([str(job)], cwd=project, SLURM_SUBMIT_DIR=str(project))
    assert result.returncode == 0, result.stderr


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
    fake.write_text('#!/bin/bash\necho "$@" > "$(dirname "$0")/called_with"\n')
    fake.chmod(0o755)
    script = (
        "set -euo pipefail\n"
        "module() { :; }\n"  # no environment modules outside the cluster
        f"source {REPO_ROOT / 'lib' / 'slurm_common.sh'}\n"
        f"source {SCRIPTS / 'get_pangenes_common.sh'}\n"
        f"GET_PANGENES_ENV={env_bin} INPUT_DIR=in N_CPUS=4 TAXA_THRESHOLD=0 VARIANT=7a\n"
        "run_get_pangenes\n"
    )
    result = _bash(["-c", script], cwd=tmp_path, TAXON_REGEX=taxon_regex)
    assert result.returncode == 0, result.stderr
    assert "get_pangenes (7a) completed" in result.stdout
    assert (env_bin / "called_with").read_text().strip() == expected_args.format(env=env_bin)


def test_run_get_pangenes_needs_the_env_setting(tmp_path):
    script = (
        f"source {REPO_ROOT / 'lib' / 'slurm_common.sh'}\n"
        f"source {SCRIPTS / 'get_pangenes_common.sh'}\n"
        "INPUT_DIR=in N_CPUS=4 TAXA_THRESHOLD=0 VARIANT=7a\n"
        "run_get_pangenes\n"
    )
    result = _bash(["-c", script], cwd=tmp_path, GET_PANGENES_ENV="", PIPELINE_CONFIG=str(tmp_path / "none.sh"))
    assert result.returncode == 1
    assert "GET_PANGENES_ENV not set" in result.stderr
