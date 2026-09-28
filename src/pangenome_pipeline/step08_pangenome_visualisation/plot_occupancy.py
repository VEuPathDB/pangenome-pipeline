#!/usr/bin/env python3
"""Pangene occupancy barplot from a get_pangenes pangene_matrix.tr.tab.

Occupancy = number of genomes a pangene (cluster) is present in. One bar per occupancy value
(1..N genomes), bar height = number of pangenes, coloured by category (counts of genomes, not fractions):
  Core       present in all N genomes                                  (--core-missing 0)
  Soft-core  not core, present in at least N-1 genomes                 (--softcore-missing 1)
  Cloud      present in at most 2 genomes, and not core or soft-core   (--cloud-max 2)
  Shell      everything in between
--core-missing / --softcore-missing set how many genomes a pangene may be missing from to still be core / soft-core.
Writes <out-prefix>.png/.svg and <out-prefix>.tsv (occupancy, category, count).

Input: rows = clusters, columns = genomes, cell > 0 = cluster present in that genome.
"""

import argparse
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pangenome_pipeline.step08_pangenome_visualisation.plot_pancore import read_presence

COLORS = {"Core": "#0072B2", "Soft-core": "#E69F00", "Shell": "#009E73", "Cloud": "#CC79A7"}
ORDER = ["Core", "Soft-core", "Shell", "Cloud"]


def categorize(occ, n, core_missing=0, softcore_missing=1, cloud_max=2):
    """Category of a pangene present in occ of the n genomes.

    Core       present in all genomes, or missing from at most core_missing of them   (default 0: all n)
    Soft-core  not core, missing from at most softcore_missing genomes                (default 1: n-1)
    Cloud      not core or soft-core, present in at most cloud_max genomes            (default 2: 1 or 2)
    Shell      everything else
    """
    if softcore_missing <= core_missing:
        raise ValueError("softcore_missing must be larger than core_missing")
    if occ >= n - core_missing:
        return "Core"
    if occ >= n - softcore_missing:
        return "Soft-core"
    if occ <= cloud_max:
        return "Cloud"
    return "Shell"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--matrix", required=True, help="pangene_matrix.tr.tab (rows = clusters, columns = genomes)")
    ap.add_argument("-o", "--out-prefix", required=True, help="output prefix (writes .png/.svg/.tsv)")
    ap.add_argument(
        "--core-missing", type=int, default=0, help="core = missing from at most this many genomes (default 0)"
    )
    ap.add_argument(
        "--softcore-missing", type=int, default=1, help="soft-core = missing from at most this many genomes (default 1)"
    )
    ap.add_argument("--cloud-max", type=int, default=2, help="cloud = present in at most this many genomes (default 2)")
    ap.add_argument("--title", default="", help="plot title (none by default)")
    args = ap.parse_args()

    genomes, pres = read_presence(args.matrix)
    n = len(genomes)
    occupancy = pres.sum(axis=1)
    counts = np.array([(occupancy == k).sum() for k in range(1, n + 1)])
    cats = [categorize(k, n, args.core_missing, args.softcore_missing, args.cloud_max) for k in range(1, n + 1)]
    print(f"{pres.shape[0]} clusters x {n} genomes", file=sys.stderr)

    with open(args.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["occupancy", "category", "count"])
        for k, (c, cat) in enumerate(zip(counts, cats), start=1):
            w.writerow([k, cat, int(c)])

    fig, ax = plt.subplots(figsize=(9, 5))
    ks = np.arange(1, n + 1)
    ax.bar(ks, counts, width=0.8, color=[COLORS[c] for c in cats])
    for k, c in zip(ks, counts):
        ax.text(k, c, f"{c}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, counts.max() * 1.08)
    ax.set_xticks(ks)
    ax.set_xlabel("Pangene occupancy", fontsize=12)
    ax.set_ylabel("Number of pangenes", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    if args.title:
        ax.set_title(args.title, loc="left", fontsize=14)
    present = [c for c in ORDER if c in cats]
    ax.legend(
        [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in present],
        present,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
    )
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
