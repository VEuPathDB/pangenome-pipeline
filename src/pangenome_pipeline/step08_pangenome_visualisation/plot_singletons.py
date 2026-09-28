#!/usr/bin/env python3
"""Singletons per genome from a get_pangenes pangene_matrix.tr.tab.

A singleton is a pangene (cluster) present in exactly one genome. One horizontal bar per genome,
bar length = number of singletons in that genome, sorted with the largest on top.
Writes <out-prefix>.png/.svg and <out-prefix>.tsv (genome, singletons).

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

SINGLETON_COLOR = "#CC79A7"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--matrix", required=True, help="pangene_matrix.tr.tab (rows = clusters, columns = genomes)")
    ap.add_argument("-o", "--out-prefix", required=True, help="output prefix (writes .png/.svg/.tsv)")
    ap.add_argument("--title", default="", help="plot title (none by default)")
    args = ap.parse_args()

    genomes, pres = read_presence(args.matrix)
    print(f"{pres.shape[0]} clusters x {len(genomes)} genomes", file=sys.stderr)
    only_one = pres.sum(axis=1) == 1
    counts = pres[only_one].sum(axis=0)
    order = np.argsort(counts)
    names = [genomes[i] for i in order]
    counts = counts[order]

    with open(args.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "singletons"])
        for name, c in zip(reversed(names), reversed(counts)):
            w.writerow([name, int(c)])

    fig, ax = plt.subplots(figsize=(7, 0.42 * len(genomes) + 1.5))
    ypos = np.arange(len(names))
    ax.barh(ypos, counts, height=0.7, color=SINGLETON_COLOR)
    for y, c in zip(ypos, counts):
        ax.text(c, y, f" {c}", ha="left", va="center", fontsize=8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names)
    ax.set_xlim(0, counts.max() * 1.1)
    ax.set_xlabel("Number of singletons", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    if args.title:
        ax.set_title(args.title, loc="left", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
