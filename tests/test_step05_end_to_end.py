"""End-to-end snapshot tests for the step 05 filter (7a) and tag (7b) scripts.

tests/data/step05/expanded_annotation holds a three-genome panel (REF, S2, S3) built so every
filter rule fires at least once:
  - geneX (no evidence) overlaps main geneB in REF, present in only 2 genomes -> dropped everywhere
  - geneE is an extra copy with no valid ORF                                  -> dropped
  - geneR is a rescued transfer at 0.30 coverage                             -> dropped
  - geneH/geneI overlap, neither has evidence, geneI has fewer good transfers -> geneI dropped
  - geneG has a small overlap on geneA's last CDS segment                    -> kept (flank rule)
  - geneF's pseudogenic_transcript child                                      -> relabelled "transcript"
  - ncA is a non-coding gene                                                  -> passed through untouched

The expected outputs in tests/data/step05/expected_* are what the scripts produced before the
refactor. Any change in them is a behaviour change: if it is intended, regenerate with
    UPDATE_SNAPSHOTS=1 nox -s pytest
and review the diff to the expected_* files before committing.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from conftest import REPO_ROOT

DATA = REPO_ROOT / "tests" / "data" / "step05"
PACKAGE = "pangenome_pipeline.step05_cross_genome_annotation_integration"
JOB = "workflow/step05_cross_genome_annotation_integration/filter_final_annotation_job.sh"

PANEL = ["--strains", "genome_order.txt", "--final-dir", "expanded_annotation"]


def _run_and_compare(script, args, expected_dir, tmp_path, pattern="*"):
    """Run a step 05 script from tests/data/step05 and compare its outputs (matching pattern) to
    the snapshot in expected_dir."""
    out_dir = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, "-m", f"{PACKAGE}.{script.removesuffix('.py')}", *args, "--out-dir", str(out_dir)],
        cwd=DATA,
        env={**os.environ, "MPLBACKEND": "Agg"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr

    if os.environ.get("UPDATE_SNAPSHOTS"):
        shutil.rmtree(expected_dir, ignore_errors=True)
        expected_dir.mkdir(parents=True)
        for f in out_dir.glob(pattern):
            shutil.copy(f, expected_dir / f.name)
        pytest.skip(f"snapshots rewritten in {expected_dir}")

    produced = sorted(p.name for p in out_dir.glob(pattern))
    assert produced == sorted(p.name for p in expected_dir.iterdir())
    for name in produced:
        assert (out_dir / name).read_text() == (expected_dir / name).read_text(), f"{name} differs from snapshot"


def test_filter_final_annotation_snapshot(tmp_path):
    _run_and_compare("filter_final_annotation.py", [*PANEL, "--id-prefix", "Tst"], DATA / "expected_filter", tmp_path)


def test_tag_final_annotation_snapshot(tmp_path):
    _run_and_compare("tag_final_annotation.py", [*PANEL, "--id-prefix", "Tst"], DATA / "expected_tag", tmp_path)


def test_plot_final_annotation_qc_tables_snapshot(tmp_path):
    # The panel's own final GFFs stand in for "each strain's original GFF" - only the tables are compared.
    args = [*PANEL, "--genome-dir", "expanded_annotation", "--genome-suffix", "_final_annotation.gff"]
    _run_and_compare("plot_final_annotation_qc.py", args, DATA / "expected_final_qc", tmp_path, pattern="*.tsv")


def test_plot_filter_qc_tables_snapshot(tmp_path):
    args = [
        "--strains",
        "genome_order.txt",
        "--original-dir",
        "expanded_annotation",
        "--filtered-dir",
        "expected_filter",
    ]
    _run_and_compare("plot_filter_qc.py", args, DATA / "expected_filter_qc", tmp_path, pattern="*.tsv")


def test_snapshot_exercises_every_filter_rule():
    """Guard the fixture itself: if it is edited so a rule no longer fires, this says so."""
    reports = "".join(p.read_text() for p in (DATA / "expected_filter").glob("*_filter_report.tsv"))
    assert "present in only 2 genome(s)" in reports
    assert "extra copy with no valid ORF" in reports
    assert "rescued transfer, Liftoff coverage 0.30" in reports
    assert "neither has evidence, geneH more supported" in reports

    s3 = (DATA / "expected_filter" / "S3_final_annotation_filtered.gff").read_text()
    assert "original_ID=geneG" in s3  # flank overlap kept
    assert "\ttranscript\t" in s3 and "pseudogenic_transcript" not in s3
    assert "ID=ncA" in (DATA / "expected_filter" / "REF_final_annotation_filtered.gff").read_text()
    # geneX never overlaps anything in S2, but the verdict is panel-wide
    assert "original_ID=geneX" not in (DATA / "expected_filter" / "S2_final_annotation_filtered.gff").read_text()
    assert "geneG\twith_main" in (DATA / "expected_tag" / "S3_intersect_pairs.tsv").read_text()


def test_fixture_inputs_untouched():
    """The scripts must never write into their --final-dir."""
    names = sorted(p.name for p in Path(DATA / "expanded_annotation").iterdir())
    assert names == ["REF_final_annotation.gff", "S2_final_annotation.gff", "S3_final_annotation.gff"]


def test_filter_job_wrapper_end_to_end(run_job, tmp_path):
    """The SLURM wrapper, submitted from a project folder as a user would, reads the build's output from
    output/step05_*/expanded_annotation/, runs the package through lib/slurm_common.sh's PYTHONPATH and
    writes the same reports as the snapshot to output/step05_*/filtered_annotation/."""
    project = tmp_path / "project"
    step_out = project / "output" / "step05_cross_genome_annotation_integration"
    shutil.copytree(DATA / "expanded_annotation", step_out / "expanded_annotation")
    shutil.copy(DATA / "genome_order.txt", project / "genome_order.txt")

    result = run_job(JOB, project)
    assert result.returncode == 0, result.stderr
    out = step_out / "filtered_annotation"
    for name in ("filter_summary.tsv", "REF_filter_report.tsv", "S2_filter_report.tsv", "S3_filter_report.tsv"):
        assert (out / name).read_text() == (DATA / "expected_filter" / name).read_text(), name


def test_check_final_annotations_reads_build_logs_from_project_logs(run_job, tmp_path):
    """A genome without its final GFF is reported with the end of its build job's log from logs/, and
    the resubmit command names the build job by its full path."""
    project = tmp_path / "project"
    final = project / "output" / "step05_cross_genome_annotation_integration" / "expanded_annotation"
    (final / "G2").mkdir(parents=True)
    (final / "G1_final_annotation.gff").write_text("##gff-version 3\n")
    (final / "G2" / "step_main.gff").write_text("##gff-version 3\n")
    (project / "genome_order.txt").write_text("G1\nG2\n")
    (project / "logs").mkdir()
    (project / "logs" / "build_final_annotation_4242.out").write_text("Array task 2: building G2\n")
    (project / "logs" / "build_final_annotation_4242.err").write_text("ERROR: AGAT failed on G2\n")

    result = run_job("workflow/step05_cross_genome_annotation_integration/check_final_annotations_job.sh", project)
    assert result.returncode == 1
    assert "OK       G1" in result.stdout
    assert "MISSING  G2  (array task 2) - stopped in round 1 (good transfers)" in result.stdout
    assert f"job log: {project}/logs/build_final_annotation_4242.err" in result.stdout
    assert "ERROR: AGAT failed on G2" in result.stdout
    assert (
        f"sbatch --array=2 {REPO_ROOT}/workflow/step05_cross_genome_annotation_integration/build_final_annotation_basic.sh"
        in result.stdout
    )
