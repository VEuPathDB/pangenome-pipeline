"""The workflow/ + docs/ + $OUTPUT_DIR layout, checked over every tracked file.

- the old numbered step folders are gone and nothing outside results/ still names them;
- every SLURM job finds the pipeline with the standard preamble and logs to logs/;
- no job writes a relative output folder into the submit folder (outputs go under $OUTPUT_DIR);
- every relative Markdown link in the pipeline's own docs points at a tracked file or folder.
"""

import re
import subprocess
from pathlib import Path

from conftest import REPO_ROOT

OLD_FOLDERS = (
    "01_Genome_Selection_QC",
    "02_TE_annotation_Clean_Up",
    "03_denovo_Genome_Annotation",
    "04_InterProScan",
    "05_Cross-genome_annotation_integration",
    "06_GET_PANGENES",
    "07_Post_GET_PANGENES_adjusting",
    "08_Pangenome_Visualisation",
    "09_Redundant_Contig_Removal",
    "Alternative_Get_PanLoci_pipeline",
    "scripts_no_slurm",
)
OLD_NAME_RE = re.compile(r"(?<![A-Za-z0-9_/])(" + "|".join(re.escape(f) for f in OLD_FOLDERS) + r")\b")

# Files allowed to name an old folder, and why.
OLD_NAMES_ALLOWED = {
    "docs/migrating_projects.md": "maps every old path to its new one",
    "scripts/sync_species_results.py": "describes the species projects' SUMMARY.md convention, which uses the old names",
    "tests/test_references.py": "KNOWN_MISSING lists never-committed files by the paths the Tgondii report uses",
    "tests/test_migrating_projects.py": "builds a project in the old layout to test the move script",
    "tests/test_workflow_layout.py": "this test",
}

PREAMBLE_RE = re.compile(
    r'_pd="\$\{PIPELINE_DIR:-\$\(cd "\$\(dirname "\$\{BASH_SOURCE\[0\]\}"\)/(?:\.\./)*\.\." && pwd\)\}"\n'
    r'\[\[ -f "\$_pd/lib/slurm_common\.sh" \]\] \|\| \{ _pd="\$\{SLURM_SUBMIT_DIR:-\$PWD\}"; '
    r'while \[\[ "\$_pd" != / && ! -f "\$_pd/lib/slurm_common\.sh" \]\]; do _pd="\$\(dirname "\$_pd"\)"; done; \}\n'
    r'source "\$_pd/lib/slurm_common\.sh" \|\| '
)

# `NAME="relative/path"` for an output-folder setting
RELATIVE_OUT_RE = re.compile(
    r'^\s*(?:local\s+)?([A-Z_]*(?:OUT|OUTDIR|OUT_DIR|OUTPUT_DIR))="(?![$/"])([^"]*)"', re.MULTILINE
)
# ...which is fine when the script first moves into its step folder (step 03, the Get-PanLoci run folder)
CD_INTO_STEP_RE = re.compile(r'^(?:.*&& )?cd "\$\{?(STEP0\d_OUT|RUN_DIR)\}?"', re.MULTILINE)
RELATIVE_OUT_ALLOWED = {
    "workflow/step05_cross_genome_annotation_integration/build_final_annotation.sh": (
        "a default only; its one caller, build_final_annotation_basic.sh, always passes -o $STEP05_OUT/..."
    ),
}

LINK_RE = re.compile(r"\]\(([^)\s#]+)(?:#[^)\s]*)?\)")


def _tracked(*patterns):
    return subprocess.run(
        ["git", "ls-files", *patterns], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()


def _jobs():
    return [f for f in _tracked("workflow/*", "species_specific/*") if f.endswith((".sh", ".slurm"))]


def test_old_step_folders_are_gone():
    assert not [f for f in OLD_FOLDERS if (REPO_ROOT / f).exists()]


def test_nothing_names_an_old_step_folder():
    hits = []
    for f in _tracked():
        if f.startswith(("results/", "tests/data/")) or f in OLD_NAMES_ALLOWED:
            continue
        path = REPO_ROOT / f
        if not path.is_file():
            continue
        try:
            text = path.read_text()
        except UnicodeDecodeError:
            continue
        hits += [
            f"{f}:{n}: {m.group(1)}" for n, line in enumerate(text.splitlines(), 1) for m in OLD_NAME_RE.finditer(line)
        ]
    assert not hits, "Old step folder names (the layout is workflow/ + docs/ now):\n  " + "\n  ".join(hits)


def test_every_slurm_job_uses_the_preamble_and_logs_to_logs():
    bad = []
    for f in _jobs():
        text = (REPO_ROOT / f).read_text()
        if not re.search(r"^#SBATCH", text, re.MULTILINE):
            continue  # a helper run with bash or sourced by a job
        if not PREAMBLE_RE.search(text):
            bad.append(f"{f}: no standard preamble (see lib/slurm_common.sh)")
        for opt in re.findall(r"^#SBATCH --(?:output|error)=(\S+)", text, re.MULTILINE):
            if not opt.startswith("logs/"):
                bad.append(f"{f}: SLURM log {opt} is not in logs/")
    assert not bad, "\n  ".join(["Jobs off the standard layout:", *bad])


def test_no_job_writes_a_relative_output_folder():
    bad = []
    for f in _jobs():
        text = (REPO_ROOT / f).read_text()
        if f in RELATIVE_OUT_ALLOWED or CD_INTO_STEP_RE.search(text):
            continue
        bad += [f'{f}: {m.group(1)}="{m.group(2)}"' for m in RELATIVE_OUT_RE.finditer(text)]
    assert not bad, "Put these under $STEPNN_OUT / $OUTPUT_DIR:\n  " + "\n  ".join(bad)


def test_markdown_links_resolve():
    tracked = set(_tracked())
    tracked_dirs = {str(p) for f in tracked for p in Path(f).parents}
    missing = []
    for f in _tracked("*.md"):
        if f.startswith("results/"):
            continue  # the species reports link their own (partly unpublished) working data
        for n, line in enumerate((REPO_ROOT / f).read_text().splitlines(), 1):
            for target in LINK_RE.findall(line):
                if "://" in target or target.startswith("mailto:"):
                    continue
                resolved = (REPO_ROOT / Path(f).parent / target).resolve()
                try:
                    rel = str(resolved.relative_to(REPO_ROOT))
                except ValueError:
                    missing.append(f"{f}:{n}: {target} (outside the repository)")
                    continue
                if rel not in tracked and rel not in tracked_dirs and rel != ".":
                    missing.append(f"{f}:{n}: {target}")
    assert not missing, "Links to files that don't exist:\n  " + "\n  ".join(missing)
