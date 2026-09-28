#!/usr/bin/env python3
"""POCS heatmap with dendrogram from the POCS.matrix.tab written by get_pangenes.

POCS = percentage of conserved sequences between two genomes. Genomes are clustered with
complete-linkage hierarchical clustering on Euclidean distances between their POCS rows.
The colour scale covers the off-diagonal values only (the diagonal, a genome against itself, is always 100
and would compress the real differences); the diagonal cells show the darkest colour and are still labelled 100.
Writes <out-prefix>.png/.svg.

Needs seaborn and pandas.
"""
import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import seaborn as sns


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--pocs", required=True, help="POCS.matrix.tab (tab-separated, genomes as rows and columns)")
    ap.add_argument("-o", "--out-prefix", required=True, help="output prefix (writes .png/.svg)")
    ap.add_argument("--cmap", default="mako_r", help="matplotlib colormap (default mako_r: darker = more similar)")
    ap.add_argument("--full-scale", action="store_true",
                    help="colour scale up to 100 (default: limited to the off-diagonal range, so the diagonal does not compress it)")
    ap.add_argument("--title", default="", help="plot title (none by default)")
    args = ap.parse_args()

    pocs = pd.read_csv(args.pocs, sep="\t", index_col=0)
    pocs.index = pocs.index.str.replace("_for_get_pangenes", "", regex=False)
    pocs.columns = pocs.columns.str.replace("_for_get_pangenes", "", regex=False)
    print(f"{pocs.shape[0]} genomes", file=sys.stderr)

    vmin = vmax = None
    if not args.full_scale:
        off = pocs.values[~np.eye(len(pocs), dtype=bool)]
        vmin, vmax = off.min(), off.max()
    g = sns.clustermap(pocs, method="complete", metric="euclidean", cmap=args.cmap, vmin=vmin, vmax=vmax,
                       figsize=(13, 14), annot=True, fmt=".0f", annot_kws={"size": 15})
    g.ax_cbar.set_xlabel("POCS")
    g.ax_heatmap.tick_params(axis="both", labelsize=15)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    if args.title:
        g.fig.suptitle(args.title, fontsize=16, x=0.02, ha="left")
    for ext in ("png", "svg"):
        g.fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out_prefix}.png/.svg", file=sys.stderr)


if __name__ == "__main__":
    main()
