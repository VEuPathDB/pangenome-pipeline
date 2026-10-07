"""No personal or cluster-specific paths in the pipeline code.

Cluster paths belong in the project's pipeline_config.sh (template: pipeline_config.example.sh,
the only code file allowed to mention them, as examples). results/ is the species reports'
record and is left alone.
"""

import re
import subprocess

from conftest import REPO_ROOT

ALLOWED = {"pipeline_config.example.sh"}
FORBIDDEN = re.compile(r"aludwig|hc-storage|/(?:home|Users|mnt|scratch|lustre|nobackup)/")


def test_no_hardcoded_paths_in_code():
    files = subprocess.run(
        ["git", "ls-files", "*.py", "*.sh", "*.slurm"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    hits = []
    for f in files:
        if f in ALLOWED or f.startswith(("results/", "tests/")):
            continue
        for lineno, line in enumerate((REPO_ROOT / f).read_text(errors="replace").splitlines(), start=1):
            if FORBIDDEN.search(line):
                hits.append(f"{f}:{lineno}: {line.strip()[:100]}")
    assert not hits, "Move these to pipeline_config.sh:\n  " + "\n  ".join(hits)
