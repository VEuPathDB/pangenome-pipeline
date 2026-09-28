"""Finding and reading GFF files."""

import os
from pathlib import Path


def find_gff(gff_dir, name, exts=(".gff3", ".gff"), case_insensitive=False):
    """Path of <gff_dir>/<name><ext> for the first ext that exists, or None.

    With case_insensitive, falls back to any *.gff* file whose stem matches name ignoring case.
    Returns a Path if gff_dir is a Path, otherwise a string, so callers keep the type they use.
    """
    as_path = isinstance(gff_dir, Path)
    for ext in exts:
        candidate = os.path.join(gff_dir, f"{name}{ext}")
        if os.path.isfile(candidate):
            return Path(candidate) if as_path else candidate
    if case_insensitive:
        lowered = name.lower()
        for f in Path(gff_dir).glob("*.gff*"):
            if f.stem.lower() == lowered:
                return f if as_path else str(f)
    return None


def gff_seqids(path):
    """Every column-1 sequence name in a GFF (comment and blank lines skipped)."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ids.add(line.split("\t", 1)[0])
    return ids
