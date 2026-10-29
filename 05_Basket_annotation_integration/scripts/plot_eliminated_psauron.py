#!/usr/bin/env python3
"""
PSAURON histograms of the genes eliminated by the step 05 quality rules, one panel per genome (4 columns, 10 bins of
0.1, the counts written on the bins from 0.7 up, the step 05 cut-offs dashed), from plot_gene_fate.py's
<label>_eliminated.tsv:
  <label>_eliminated_psauron.png        every eliminated gene with a PSAURON score
  <label>_eliminated_psauron_only.png   the genes whose only sign is PSAURON (the "PSAURON only" bar of
                                        <label>_eliminated_reasons.png)

Usage:
    python3 plot_eliminated_psauron.py --eliminated basket_annotation/plots/Toxo_prefer_panel_model_eliminated.tsv \
        --label Toxo_prefer_panel_model --out-dir basket_annotation/plots
"""
import argparse
import csv
import os


def hist(rows, genomes, a, what, name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    nc = 4
    nr = -(-len(genomes) // nc)
    fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 2.3 * nr), sharex=True, sharey=True, squeeze=False)
    bins = np.linspace(0, 1, 11)   # 10 bins of 0.1
    n_all = 0
    for axx, g in zip(axs.flat, genomes):
        ps = [float(r["psauron"]) for r in rows if r["genome"] == g and r["psauron"] not in ("", "-")]
        n_all += len(ps)
        n, _, _ = axx.hist(ps, bins=bins, color="#2a6f97", edgecolor="white", lw=0.5)
        for lo, k in zip(bins[:-1], n):   # counts on the bins from 0.7 up (the borderline ones)
            if lo >= 0.7 - 1e-9:
                axx.text(lo + 0.05, k, str(int(k)), ha="center", va="bottom", fontsize=7.5)
        for t in (a.psauron_min, a.psauron_min_in_ref):
            axx.axvline(t, color="#c0392b", ls="--", lw=0.9)
        axx.set_title(f"{g} ({len(ps)})", fontsize=9, loc="left")
        axx.tick_params(labelbottom=True)   # the PSAURON scale on every panel, not only the bottom row
        axx.spines[["top", "right"]].set_visible(False)
    for axx in list(axs.flat)[len(genomes):]:
        axx.set_visible(False)
    fig.supxlabel(f"PSAURON (dashed: step 05 cut-offs {a.psauron_min} and {a.psauron_min_in_ref})", fontsize=10)
    fig.supylabel("genes eliminated", fontsize=10)
    fig.suptitle(f"{a.label}: PSAURON of {what}, by genome ({n_all:,} genes)", x=0.01, ha="left", fontsize=12)
    fig.tight_layout()
    out = os.path.join(a.out_dir, name)
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eliminated", required=True, help="plot_gene_fate.py's <label>_eliminated.tsv")
    ap.add_argument("--label", default="panel")
    ap.add_argument("--psauron-min", type=float, default=0.7)
    ap.add_argument("--psauron-min-in-ref", type=float, default=0.9)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rows = list(csv.DictReader(open(a.eliminated), delimiter="\t"))
    genomes = list(dict.fromkeys(r["genome"] for r in rows))
    hist(rows, genomes, a, "the genes eliminated by the step 05 quality rules", f"{a.label}_eliminated_psauron.png")
    hist([r for r in rows if r["signs"] == "low_psauron"], genomes, a, "the genes eliminated by PSAURON only",
         f"{a.label}_eliminated_psauron_only.png")


if __name__ == "__main__":
    main()
