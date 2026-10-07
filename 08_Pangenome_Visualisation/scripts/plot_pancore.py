#!/usr/bin/env python3
"""Pan/core genome curves from a get_pangenes pangene_matrix.tr.tab.

Python port of plot_pancore_matrix.pl (GET-HOMOLOGUES, Contreras-Moreira & Vinuesa).
Genomes are added one by one in random order (many permutations). For each step:
  pan  = clusters present in at least one genome added so far
  core = clusters present in all genomes added so far
Fits: pan (Tettelin 2005), core Tettelin (2005) and core Willenbrock (2007).

Input: rows = clusters, columns = genomes, cell > 0 = cluster present in that genome.
Writes <out-prefix>.png/.svg (pan and core in one panel: all permutation points + fitted curves) and <out-prefix>.log.
"""

import argparse
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

TAU_START, TAU_STEP = 2.8, 0.2

PAN_COLOR, CORE_COLOR = "#E69F00", "#0072B2"


def read_presence(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    keep = [i for i, g in enumerate(rows[0]) if i > 0 and g.strip()]
    genomes = [rows[0][i].removesuffix("_for_get_pangenes") for i in keep]
    pres = []
    for r in rows[1:]:
        vals = [(i < len(r) and r[i] not in ("", "NA", "0", "-")) for i in keep]
        if any(vals):
            pres.append(vals)
    return genomes, np.array(pres, dtype=bool)


def permutations(pres, n_perm, seed):
    rng = np.random.default_rng(seed)
    n_genomes = pres.shape[1]
    xs, pan, core = [], [], []
    for _ in range(n_perm):
        order = rng.permutation(n_genomes)
        sub = pres[:, order].cumsum(axis=1)
        pan.append((sub > 0).sum(axis=0))
        core.append((sub == np.arange(1, n_genomes + 1)).sum(axis=0))
        xs.append(np.arange(1, n_genomes + 1))
    return np.concatenate(xs), np.concatenate(pan), np.concatenate(core)


def core_tettelin(g, omega, k, tau):
    return omega + k * np.exp(-g / tau)


def core_willenbrock(g, omega, k, tau):
    return omega + k * np.exp(-np.sqrt(np.sqrt(g)) / tau)


def pan_tettelin(g, omega, theta, k, tau):
    return omega + theta * (g - 1) + k * np.exp(-2 / tau) * ((1 - np.exp(-(g - 1) / tau)) / (1 - np.exp(-1 / tau)))


def fit(func, x, y, p0_fn):
    """Try decreasing tau starting values until nls-style fitting converges, as the perl does."""
    tau = TAU_START + TAU_STEP
    while tau > 0:
        try:
            popt, _ = curve_fit(func, x, y, p0=p0_fn(tau), maxfev=20000)
            resid = y - func(x, *popt)
            ser = float(np.sqrt((resid**2).sum() / (len(y) - len(popt))))
            return popt, ser
        except (RuntimeError, ValueError):
            tau -= TAU_STEP
    return None, None


def summarize(x, y, n_genomes):
    ks = np.arange(1, n_genomes + 1)
    rows = np.array([[y[x == g].mean(), y[x == g].min(), y[x == g].max()] for g in ks])
    return ks, rows[:, 0], rows[:, 1], rows[:, 2]


def plot_curves(ax, x, series, fits, n_genomes, title, style="points"):
    """series = [(y, color, label), ...] drawn as dots or mean+band; fits = [(func, popt, color, ls, label), ...]."""
    xs = np.linspace(1, n_genomes, 200)
    for y, color, label in series:
        if style == "points":
            ax.scatter(x, y, s=10, color=color, alpha=0.35, linewidths=0)
        else:
            ks, mean, lo, hi = summarize(x, y, n_genomes)
            ax.fill_between(ks, lo, hi, color=color, alpha=0.15, linewidth=0)
            ax.plot(ks, mean, color=color, lw=2.2, marker="o", ms=4.5, label=label)
    for func, popt, color, ls, label in fits:
        if popt is not None:
            ax.plot(
                xs,
                func(xs, *popt),
                color=color,
                lw=2.2 if style == "points" else 1.6,
                ls=ls if style == "points" else "--",
                label=label,
            )
    ax.set_xlabel("Number of genomes sampled", fontsize=12)
    ax.set_ylabel("Number of pangenes", fontsize=12)
    ax.set_xticks(range(1, n_genomes + 1))
    ax.spines[["top", "right"]].set_visible(False)
    if title:
        ax.set_title(title, loc="left", fontsize=14, pad=30)
    handles, labels = ax.get_legend_handles_labels()
    for h in handles:
        if hasattr(h, "set_alpha"):
            h.set_alpha(1)
    ax.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=3, frameon=False, fontsize=9)


def plot_panel(ax, x, y, ylabel, curves, n_genomes):
    ax.scatter(x, y, s=6, color="0.35", alpha=0.35, linewidths=0)
    xs = np.linspace(1, n_genomes, 200)
    for func, popt, color, label in curves:
        if popt is not None:
            ax.plot(xs, func(xs, *popt), color=color, lw=2, label=label)
    ax.set_xlabel("genomes (g)")
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(1, n_genomes + 1))
    ax.spines[["top", "right"]].set_visible(False)
    if any(p is not None for _, p, _, _ in curves):
        ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=2, frameon=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--matrix", required=True, help="pangene_matrix.tr.tab (rows = clusters, columns = genomes)")
    ap.add_argument("-o", "--out-prefix", required=True, help="output prefix (writes .png/.svg/.log)")
    ap.add_argument("--permutations", type=int, default=100, help="random genome orders to sample (default 100)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument(
        "--title",
        default="",
        help='plot title (combined layout: default "Pangenome and core genome size"; separate layout: overall title above both panels, none by default)',
    )
    ap.add_argument(
        "--layout",
        choices=["combined", "separate"],
        default="combined",
        help="combined = pan and core in one panel; separate = pan and core in two panels, each with its own y-axis (ignored for --style panels)",
    )
    ap.add_argument(
        "--style",
        choices=["points", "mean", "panels"],
        default="points",
        help="points = every permutation as a dot; mean = mean line + min-max band; panels = original .pl layout, pan and core in separate panels with grey dots (default points)",
    )
    args = ap.parse_args()

    genomes, pres = read_presence(args.matrix)
    n = len(genomes)
    print(f"{pres.shape[0]} clusters x {n} genomes", file=sys.stderr)
    x, pan, core = permutations(pres, args.permutations, args.seed)
    xf = x.astype(float)

    pan_p, pan_ser = fit(pan_tettelin, xf, pan.astype(float), lambda tau: [pan.min(), 33.0, 476.0, tau])
    core_t, core_t_ser = fit(
        core_tettelin, xf, core.astype(float), lambda tau: [core.min(), core.max() - core.min(), tau]
    )
    core_w, core_w_ser = fit(
        core_willenbrock, xf, core.astype(float), lambda tau: [core.min(), core.max() - core.min(), tau]
    )

    with open(args.out_prefix + ".log", "w") as log:
        log.write(f"# {pres.shape[0]} clusters, {n} genomes, {args.permutations} permutations, seed {args.seed}\n")
        for name, popt, ser in (
            ("pan (Tettelin)", pan_p, pan_ser),
            ("core Tettelin", core_t, core_t_ser),
            ("core Willenbrock", core_w, core_w_ser),
        ):
            if popt is None:
                log.write(f"# {name}: fit failed to converge\n")
            else:
                log.write(f"# {name}: params = {', '.join(f'{v:.3f}' for v in popt)}; residual SE = {ser:.2f}\n")
        log.write("# genomes\tmean_pan\tmean_core\n")
        for g in range(1, n + 1):
            m = x == g
            log.write(f"{g}\t{pan[m].mean():.1f}\t{core[m].mean():.1f}\n")

    if args.style == "panels":
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
        plot_panel(axes[0], x, pan, "pan genome size (clusters)", [(pan_tettelin, pan_p, "tab:red", "Tettelin fit")], n)
        plot_panel(
            axes[1],
            x,
            core,
            "core genome size (clusters)",
            [
                (core_tettelin, core_t, "tab:red", "Tettelin fit"),
                (core_willenbrock, core_w, "tab:blue", "Willenbrock fit"),
            ],
            n,
        )
    else:
        pan_fits = [(pan_tettelin, pan_p, PAN_COLOR, "-", "Pan fit (Tettelin)")]
        core_fits = [
            (core_tettelin, core_t, CORE_COLOR, "-", "Core fit (Tettelin)"),
            (core_willenbrock, core_w, "#56B4E9", "--", "Core fit (Willenbrock)"),
        ]
        pan_series = (pan, PAN_COLOR, "Pangenome")
        core_series = (core, CORE_COLOR, "Core genome")
        if args.layout == "separate":
            fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
            plot_curves(axes[0], x, [pan_series], pan_fits, n, "Pangenome size", args.style)
            plot_curves(axes[1], x, [core_series], core_fits, n, "Core genome size", args.style)
            if args.title:
                fig.suptitle(args.title, fontsize=15, x=0.01, ha="left")
        else:
            fig, ax = plt.subplots(figsize=(8, 5.5))
            plot_curves(
                ax,
                x,
                [pan_series, core_series],
                pan_fits + core_fits,
                n,
                args.title or "Pangenome and core genome size",
                args.style,
            )
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out_prefix}.png/.svg/.log", file=sys.stderr)


if __name__ == "__main__":
    main()
