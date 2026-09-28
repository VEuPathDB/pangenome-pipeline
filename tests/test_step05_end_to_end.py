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
SCRIPTS = REPO_ROOT / "05_Cross-genome_annotation_integration" / "scripts"


def _run_and_compare(script, expected_dir, tmp_path):
    out_dir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / script),
            "--strains",
            "genome_order.txt",
            "--final-dir",
            "expanded_annotation",
            "--out-dir",
            str(out_dir),
            "--id-prefix",
            "Tst",
        ],
        cwd=DATA,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr

    if os.environ.get("UPDATE_SNAPSHOTS"):
        shutil.rmtree(expected_dir, ignore_errors=True)
        shutil.copytree(out_dir, expected_dir)
        pytest.skip(f"snapshots rewritten in {expected_dir}")

    produced = sorted(p.name for p in out_dir.iterdir())
    assert produced == sorted(p.name for p in expected_dir.iterdir())
    for name in produced:
        assert (out_dir / name).read_text() == (expected_dir / name).read_text(), f"{name} differs from snapshot"


def test_filter_final_annotation_snapshot(tmp_path):
    _run_and_compare("filter_final_annotation.py", DATA / "expected_filter", tmp_path)


def test_tag_final_annotation_snapshot(tmp_path):
    _run_and_compare("tag_final_annotation.py", DATA / "expected_tag", tmp_path)


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
