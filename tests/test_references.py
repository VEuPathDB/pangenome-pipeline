"""Every script named in a README, wrapper or docstring must exist in the repo.

This is the safety net for deleting, merging or renaming scripts: any README link, `sbatch ...`
example, "run after X.py" comment or wrapper path left pointing at the old name fails here.

A reference with a directory part (`scripts/foo.py`, `05_.../scripts/foo.py`) must resolve
relative to the referencing file, the repo root, or match the tail of a tracked path. A bare
name (`foo.py`) must be the basename of some tracked file.
"""

import re
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Third-party tools the docs mention by filename - not part of this repo - and pipeline_config.sh,
# which each user creates from pipeline_config.example.sh and is never committed.
EXTERNAL = {"MCHelper.py", "famdb.py", "interproscan.sh", "pipeline_config.sh"}

# Referenced but deliberately absent. Each needs a reason, and test_known_missing_are_still_missing
# makes sure an entry is removed once the file shows up.
KNOWN_MISSING = {
    "01_Genome_Selection_QC/scripts/contig_redundancy_final.slurm.sh": (
        "Tgondii-only script kept in the species project, not committed (refactor plan A4)"
    ),
    "contig_redundancy_final.slurm.sh": "same as above, mentioned by name in step 09's make_removal_list.py",
    "09_Redundant_Contig_Removal/scripts/remove_redundant_contig_genes.slurm.sh": (
        "Tgondii-only step 09 wrapper, not committed (refactor plan A4)"
    ),
    "09_Redundant_Contig_Removal/scripts/plot_nonredundant_job.sh": (
        "Tgondii-only step 09 wrapper, not committed (refactor plan A4)"
    ),
    "plot_nonredundant_job.sh": "same as above, mentioned by name in the Tgondii report",
    "01_Genome_Selection_QC/scripts/summarize_dropped_genes.py": (
        "Tgondii-only helper called by species_specific/Tgondii/.../drop_short_contigs.slurm.sh, not committed"
    ),
}

SCANNED_PATTERNS = ("*.md", "*.sh", "*.slurm", "*.py")
SCRIPT_EXT = r"\.(?:py|sh|slurm)"
TOKEN_RE = re.compile(rf"[A-Za-z0-9_./\-]*[A-Za-z0-9_\-]{SCRIPT_EXT}(?![A-Za-z0-9_])")
SEGMENT_IS_SCRIPT_RE = re.compile(rf"{SCRIPT_EXT}$")

# A token right after one of these is a template or glob, e.g. $NAME.sh, <genome>.py, *.slurm.sh
PLACEHOLDER_PREFIXES = ("$", "{", "<", "*", "%")
# A leading "/" right after one of these belongs to a variable, e.g. ${SLURM_SUBMIT_DIR}/05_.../x.py
VARIABLE_PREFIXES = ("}", ")")


def _git_ls_files(*patterns):
    return subprocess.run(
        ["git", "ls-files", *patterns], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()


TRACKED = set(_git_ls_files())
TRACKED_NAMES = {Path(f).name for f in TRACKED}


def _split_joined(token):
    """'a.py/b.py' (two names joined by a slash in prose) -> ['a.py', 'b.py']."""
    parts, current = [], []
    for segment in token.split("/"):
        current.append(segment)
        if SEGMENT_IS_SCRIPT_RE.search(segment):
            parts.append("/".join(current))
            current = []
    return parts or [token]


def extract_references(text):
    """Script paths/names referenced in text, normalised, placeholders and absolute paths skipped."""
    refs = []
    for match in TOKEN_RE.finditer(text):
        token = match.group(0)
        before = text[match.start() - 1] if match.start() else ""
        if before == "$" and "/" in token:
            token = token.split("/", 1)[1]  # $PROJECT_DIR/05_.../x.py -> check the part after the variable
        elif before in PLACEHOLDER_PREFIXES:
            continue
        if token.startswith("/"):
            if before not in VARIABLE_PREFIXES:
                continue  # absolute path outside the repo - Phase C's hardcoded-path check covers these
            token = token.lstrip("/")
        token = token.removeprefix("./")
        if token.startswith("..."):  # "…/foo.sh" shorthand in prose
            token = token.lstrip(".").lstrip("/")
        if not Path(token).stem or token.startswith("."):  # a bare extension like ".slurm.sh"
            continue
        refs.extend(_split_joined(token))
    return refs


def resolves(ref, referencing_file):
    if Path(ref).name in EXTERNAL:
        return True
    if "/" not in ref:
        return ref in TRACKED_NAMES
    for base in (REPO_ROOT / Path(referencing_file).parent, REPO_ROOT):
        try:
            if str((base / ref).resolve().relative_to(REPO_ROOT)) in TRACKED:
                return True
        except ValueError:
            pass
    return any(f == ref or f.endswith("/" + ref) for f in TRACKED)


def find_missing():
    missing = []
    for path in _git_ls_files(*SCANNED_PATTERNS):
        if path.startswith("tests/"):
            continue  # the tests name made-up scripts on purpose
        if path.startswith("results/") and not path.endswith(".md"):
            continue  # frozen copies of old script versions; the reports themselves are still checked
        text = (REPO_ROOT / path).read_text(errors="replace")
        for lineno, line in enumerate(text.splitlines(), start=1):
            for ref in extract_references(line):
                if ref not in KNOWN_MISSING and not resolves(ref, path):
                    missing.append(f"{path}:{lineno}: {ref}")
    return missing


def test_every_referenced_script_exists():
    missing = find_missing()
    assert not missing, "References to scripts that don't exist:\n  " + "\n  ".join(missing)


def test_known_missing_are_still_missing():
    now_present = [ref for ref in KNOWN_MISSING if resolves(ref, "README.md")]
    assert not now_present, f"These now exist - remove them from KNOWN_MISSING: {now_present}"


#
# The scanner itself
#


def test_extract_references_normalises_tokens():
    line = (
        "sbatch ./05_x/scripts/run.slurm.sh; see ${SLURM_SUBMIT_DIR}/01_y/scripts/a.py, $PROJECT_DIR/02_z/e.py, "
        "then .../b.sh and c.py/d.py; not $NAME.sh, <genome>.py, *.slurm.sh, /opt/tool.py or .slurm.sh"
    )
    assert extract_references(line) == [
        "05_x/scripts/run.slurm.sh",
        "01_y/scripts/a.py",
        "02_z/e.py",
        "b.sh",
        "c.py",
        "d.py",
    ]


def test_resolves_relative_to_file_root_and_tail():
    assert resolves(
        "../workflow/step01_genome_selection_qc/drop_contigs_slurm.sh", "docs/step01_genome_selection_qc.md"
    )  # relative to the doc
    assert resolves("workflow/step01_genome_selection_qc/drop_contigs_slurm.sh", "README.md")  # relative to the root
    assert resolves("slurm/config.sh", "somewhere/else.md")  # tail of a tracked path
    assert resolves("drop_contigs.py", "anything")  # bare name
    assert not resolves("scripts/no_such_script.py", "README.md")
    assert not resolves("no_such_script.py", "README.md")
