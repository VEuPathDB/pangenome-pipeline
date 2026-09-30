#!/usr/bin/env python3
"""Soft-core pangenes missing per genome, from a get_pangenes pangene_matrix.tr.tab.

Soft-core pangenes are present in all but --softcore-missing genomes (same rule as
plot_occupancy.py). At the default cuts (core-missing 0, softcore-missing 1) each soft-core
pangene is missing from exactly one genome, so summing that across all soft-core pangenes, per
genome, gives a direct count: how many soft-core pangenes is THIS genome missing.

Soft-core genes are expected to be shared by nearly every genome in the panel, so no single
genome should carry a much larger share of that "missing" total than the others - if one does,
that is a genuine content-QC signal (not just noise), worth checking against that genome's own
assembly (dropped or short contigs, low-coverage regions) before trusting its occupancy calls.

Writes <out-prefix>.png/.svg (one bar per genome, worst first, with the panel mean as a dashed
reference line) and <out-prefix>.tsv (genome, softcore_missing, softcore_total, pct_missing).

Input: rows = clusters, columns = genomes, cell > 0 = cluster present in that genome.
"""
import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_occupancy import categorize
from plot_pancore import read_presence


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--matrix", required=True, help="pangene_matrix.tr.tab (rows = clusters, columns = genomes)")
    ap.add_argument("-o", "--out-prefix", required=True, help="output prefix (writes .png/.svg/.tsv)")
    ap.add_argument("--core-missing", type=int, default=0, help="core = missing from at most this many genomes (default 0)")
    ap.add_argument("--softcore-missing", type=int, default=1, help="soft-core = missing from at most this many genomes (default 1)")
    ap.add_argument("--cloud-max", type=int, default=2, help="cloud = present in at most this many genomes (default 2, unused here but kept for the same category cut as the other plots)")
    ap.add_argument("--title", default="", help="plot title (none by default)")
    args = ap.parse_args()

    genomes, pres = read_presence(args.matrix)
    n = len(genomes)
    occupancy = pres.sum(axis=1)
    cats = np.array([categorize(int(k), n, args.core_missing, args.softcore_missing, args.cloud_max) for k in occupancy])
    softcore_mask = cats == "Soft-core"
    n_softcore = int(softcore_mask.sum())
    print(f"{pres.shape[0]} clusters x {n} genomes, {n_softcore} soft-core", file=sys.stderr)

    missing_counts = (~pres[softcore_mask]).sum(axis=0) if n_softcore else np.zeros(n, dtype=int)
    order = np.argsort(-missing_counts)  # worst (most missing) first
    mean = float(missing_counts.mean()) if n_softcore else 0.0

    with open(args.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "softcore_missing", "softcore_total", "pct_missing"])
        for i in order:
            pct = 100 * missing_counts[i] / n_softcore if n_softcore else 0.0
            w.writerow([genomes[i], int(missing_counts[i]), n_softcore, f"{pct:.2f}"])

    fig, ax = plt.subplots(figsize=(max(9, n * 0.5), 5))
    xs = np.arange(n)
    ax.bar(xs, missing_counts[order], width=0.75, color="#E69F00")
    for x, v in zip(xs, missing_counts[order]):
        ax.text(x, v, f"{int(v)}", ha="center", va="bottom", fontsize=8)
    ax.axhline(mean, color="#555555", linestyle="--", linewidth=1, label=f"mean = {mean:.1f}")
    ax.set_xticks(xs)
    ax.set_xticklabels([genomes[i] for i in order], rotation=45, ha="right")
    ax.set_ylim(0, max(1, missing_counts.max()) * 1.15)
    ax.set_ylabel(f"Soft-core pangenes missing (of {n_softcore})", fontsize=12)
    ax.set_xlabel("Genome", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    if args.title:
        ax.set_title(args.title, loc="left", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
