"""Reading the small list/config files the pipeline passes between steps."""

import os
from pathlib import Path


def read_strains(value):
    """Genome names, in order: from a file with one per line ('#' lines and blank lines skipped),
    or from a comma-separated string if `value` isn't a file. The first one is the reference."""
    if os.path.isfile(value):
        with open(value) as fh:
            return [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    return [s.strip() for s in value.split(",") if s.strip()]


def read_id_list(value):
    """IDs from a file (first whitespace-separated column per line; '#' lines and blank lines
    skipped), or from a comma-separated string if `value` isn't a file."""
    if os.path.isfile(value):
        with open(value) as fh:
            return [line.split()[0] for line in fh if line.strip() and not line.startswith("#")]
    return [g.strip() for g in value.split(",") if g.strip()]


def load_extra_genes(path) -> set:
    """(genome_label, gene_id) pairs from a 2-column TSV (no header; '#' lines and blank lines skipped)."""
    pairs = set()
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genome_label, gene_id = line.split("\t")
            pairs.add((genome_label, gene_id))
    return pairs


def read_paf(path, min_aln):
    """minimap2 PAF alignments at least min_aln bp long, one dict per line (0-based, half-open)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            if int(f[10]) < min_aln:
                continue
            rows.append(
                {
                    "q": f[0],
                    "qlen": int(f[1]),
                    "qs": int(f[2]),
                    "qe": int(f[3]),
                    "strand": f[4],
                    "t": f[5],
                    "tlen": int(f[6]),
                    "ts": int(f[7]),
                    "te": int(f[8]),
                    "nmatch": int(f[9]),
                    "alen": int(f[10]),
                }
            )
    return rows
