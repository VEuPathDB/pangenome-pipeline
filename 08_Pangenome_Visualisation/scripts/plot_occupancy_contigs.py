#!/usr/bin/env python3
"""Size of the contigs that carry each pangene, by pangene occupancy.

Only pangenes with exactly one gene in every genome they are present in are used. Each pangene gets ONE
value: the median, over its genes, of the length of the contig the gene lies on. One boxplot per
occupancy, coloured by category (Core / Soft-core / Shell / Cloud, same definition as plot_occupancy.py).
The Tukey HSD p-value of Core versus Cloud is computed on log10 of the contig length and written on the
line joining the two; genes of one pangene come from different genomes and are not independent samples,
so it is optimistic. The .tsv also gives, per category, the share of pangenes whose median contig is
shorter than each of --short-kb thresholds.

Inputs:
  --genes-matrix   pangene_matrix_genes.tr.tab from get_pangenes (cells = gene IDs, "-" = absent)
  --gff-dir        GFFs whose gene IDs match the matrix, named <genome><--gff-suffix>
  --contig-sizes   tab-separated table with columns genome_label, contig, length
                   (step 01: contig_size_profile/contig_sizes.tsv)
Writes <out-prefix>.png/.svg and <out-prefix>.tsv (one row per pangene used).

Needs numpy, scipy, matplotlib.
"""

import argparse
import csv
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter, NullFormatter
from matplotlib.transforms import blended_transform_factory
from plot_occupancy import COLORS, ORDER, categorize
from scipy.stats import tukey_hsd

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")


def gene_contigs(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9 and f[2] == "gene":
                m = ID_RE.search(f[8])
                if m:
                    out[m.group(1)] = f[0]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_for_get_pangenes.gff")
    ap.add_argument("--contig-sizes", required=True)
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument(
        "--short-kb",
        default="50,100,500",
        help="contig-length thresholds (kb) for the share of pangenes on short contigs",
    )
    ap.add_argument("--yscale", choices=["log", "linear"], default="log")
    ap.add_argument(
        "--style",
        choices=["boxplot", "classes"],
        default="boxplot",
        help="boxplot = median contig length per occupancy; classes = share of the pangenes of each category "
        "in each contig-size class (stacked bars)",
    )
    ap.add_argument(
        "--core-missing", type=int, default=0, help="core = missing from at most this many genomes (default 0)"
    )
    ap.add_argument(
        "--softcore-missing", type=int, default=1, help="soft-core = missing from at most this many genomes (default 1)"
    )
    ap.add_argument("--cloud-max", type=int, default=2, help="cloud = present in at most this many genomes (default 2)")
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    sizes = {}
    with open(args.contig_sizes, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sizes[(r["genome_label"], r["contig"])] = int(r["length"])

    with open(args.genes_matrix, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    cols = [(i, g.replace("_for_get_pangenes", "")) for i, g in enumerate(rows[0]) if i > 0 and g.strip()]
    n_genomes = len(cols)

    contig_of = {}
    for _, genome in cols:
        path = os.path.join(args.gff_dir, genome + args.gff_suffix)
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")
        contig_of[genome] = gene_contigs(path)

    used, multi, missing = [], 0, 0
    for r in rows[1:]:
        if r[0].startswith("chr"):
            continue
        members = []
        for i, genome in cols:
            cell = r[i] if i < len(r) else ""
            if cell not in ("", "-", "NA"):
                members.append((genome, cell.split(",")))
        if not members:
            continue
        if any(len(ids) > 1 for _, ids in members):
            multi += 1
            continue
        lengths = []
        for genome, ids in members:
            c = contig_of[genome].get(ids[0])
            L = sizes.get((genome, c)) if c else None
            if L is None:
                lengths = None
                break
            lengths.append(L)
        if lengths is None:
            missing += 1
            continue
        occ = len(members)
        used.append(
            (
                r[0],
                occ,
                categorize(occ, n_genomes, args.core_missing, args.softcore_missing, args.cloud_max),
                float(np.median(lengths)),
            )
        )
    print(
        f"{len(used)} pangenes used; {multi} left out (more than one gene in a genome)"
        + (f"; {missing} skipped (contig length not found)" if missing else ""),
        file=sys.stderr,
    )
    if not used:
        sys.exit("ERROR: no pangene left; check the GFF gene IDs, the contig names and --contig-sizes")

    with open(args.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pangene", "occupancy", "category", "median_contig_bp"])
        for u in used:
            w.writerow([u[0], u[1], u[2], f"{u[3]:.0f}"])

    occ_all = np.array([u[1] for u in used])
    cat_all = np.array([u[2] for u in used])
    val = np.array([u[3] for u in used])

    print("\nShare of pangenes whose median contig is shorter than:")
    thr = [int(x) for x in args.short_kb.split(",")]
    print(f"  {'category':10s}{'pangenes':>9s}{'median contig (kb)':>20s}" + "".join(f"{f'<{t} kb':>10s}" for t in thr))
    for cat in ORDER:
        v = val[cat_all == cat]
        if len(v):
            print(
                f"  {cat:10s}{len(v):9d}{np.median(v) / 1000:20.0f}"
                + "".join(f"{100 * (v < t * 1000).mean():9.1f}%" for t in thr)
            )

    if args.style == "classes":
        classes = [
            (0, 20e3, "10-20 kb"),
            (20e3, 50e3, "20-50 kb"),
            (50e3, 200e3, "50-200 kb"),
            (200e3, 500e3, "200-500 kb"),
            (500e3, 1e6, "0.5-1 Mb"),
            (1e6, float("inf"), "> 1 Mb"),
        ]
        shades = ["#7f2704", "#a63603", "#d94801", "#fd8d3c", "#fdd0a2", "#e8e8e8"]
        cats = [
            c for c in reversed(ORDER) if (cat_all == c).any()
        ]  # Cloud on the left, Core on the right, like the occupancy plot
        fig, ax = plt.subplots(figsize=(7.5, 5.2))
        bottom = np.zeros(len(cats))
        for (lo, hi, name), col in zip(classes, shades):
            share = np.array([100 * ((val[cat_all == c] >= lo) & (val[cat_all == c] < hi)).mean() for c in cats])
            ax.bar(range(len(cats)), share, bottom=bottom, color=col, width=0.7, label=name, edgecolor="white", lw=0.5)
            for x, (sh, b) in enumerate(zip(share, bottom)):
                if sh >= 3:
                    ax.text(
                        x,
                        b + sh / 2,
                        f"{sh:.0f}%",
                        ha="center",
                        va="center",
                        fontsize=8,
                        color="white" if col in shades[:3] else "black",
                    )
            bottom += share
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels([f"{c}\n(n={int((cat_all == c).sum())})" for c in cats])
        ax.set_ylabel("Percentage of pangenes", fontsize=12)
        ax.set_ylim(0, 100)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(title="Median contig length", loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
        if args.title:
            fig.suptitle(args.title, x=0.01, ha="left", fontsize=14)
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
        print(f"Wrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)
        return

    groups = [np.log10(val[cat_all == c]) for c in ORDER if (cat_all == c).sum() > 1]
    names = [c for c in ORDER if (cat_all == c).sum() > 1]
    p = None
    if "Core" in names and "Cloud" in names:
        p = float(tukey_hsd(*groups).pvalue[names.index("Core"), names.index("Cloud")])

    ks = [k for k in range(1, n_genomes + 1) if (occ_all == k).any()]
    cat_of = {k: cat_all[occ_all == k][0] for k in ks}
    xs = np.arange(len(ks))
    fig, ax = plt.subplots(figsize=(9, 4.8))
    bp = ax.boxplot(
        [val[occ_all == k] for k in ks],
        positions=xs,
        widths=0.7,
        patch_artist=True,
        medianprops=dict(color="black", lw=1.5),
        flierprops=dict(marker="o", ms=2.5, alpha=0.4, mec="none", mfc="0.3"),
    )
    for patch, k in zip(bp["boxes"], ks):
        patch.set_facecolor(COLORS[cat_of[k]])
        patch.set_alpha(0.9)
    plain = FuncFormatter(lambda v, _: f"{v:,.0f}")
    if args.yscale == "log":
        ax.set_yscale("log")
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(top=val.max() * 4)
    else:
        ax.set_ylim(0, val.max() * 1.25)
    ax.yaxis.set_major_formatter(plain)
    ax.set_ylabel("Contig length (bp),\nmedian per pangene", fontsize=12)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{k}\n(n={int((occ_all == k).sum())})" for k in ks], fontsize=8)
    ax.set_xlabel("Pangene occupancy", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)

    cloud_x = [x for x, k in zip(xs, ks) if cat_of[k] == "Cloud"]
    core_x = [x for x, k in zip(xs, ks) if cat_of[k] == "Core"]
    if p is not None and cloud_x and core_x:
        t = blended_transform_factory(ax.transData, ax.transAxes)
        x0, x1 = float(np.mean(cloud_x)), float(np.mean(core_x))
        ax.plot([x0, x0, x1, x1], [0.90, 0.93, 0.93, 0.90], transform=t, color="black", lw=1, clip_on=False)
        ax.text(
            (x0 + x1) / 2,
            0.94,
            "p < 0.001" if p < 0.001 else f"p = {p:.3f}",
            transform=t,
            ha="center",
            va="bottom",
            fontsize=9,
        )
    present = [c for c in ORDER if c in cat_of.values()]
    ax.legend(
        [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in present],
        present,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.12),
        ncol=len(present),
        frameon=False,
    )
    if args.title:
        fig.suptitle(args.title, x=0.01, ha="left", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"\nWrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
