#!/usr/bin/env python3
"""
Pilot: does the number of orthologs (VEuPathDB OrthoMCL) of a pangene go with its quality?

Per pangene (explorer data with orthogroups, from add_orthogroups.py):
  orthologs      the "Ortholog count" of its reference gene (largest one if it has several; none if no reference gene)
  psauron        PSAURON of its panel gene model (step 05)
  share_invalid  strains whose genes all have an invalid ORF / judged strains
  low_confidence share_invalid > 0.5 and not a reference pseudogene (same rule as the explorer)

Outputs (--out-dir): <label>_orthologs_vs_quality.png, <label>_orthologs_invalid_psauron_grid_<N>plus_strains.png (the three together, pangenes in N+ strains:
orthologs x share of invalid strains, each cell = pangenes, median PSAURON and % with PSAURON < 0.5) and
<label>_orthologs_psauron_invalid_scatter_<N>plus_strains.png (one point per pangene: orthologs x PSAURON,
coloured by the share of invalid strains, with the Spearman correlations of each pair and partial ones) and
orthologs_vs_quality.tsv (one row per pangene)

Usage:
    python3 plot_orthologs_vs_quality.py --explorer-dir pangenome_explorer/preliminary_orthogroups \\
        --og-table GenesByOrthologPattern_Summary.txt --label Toxo --out-dir pangenome_explorer/preliminary_orthogroups
"""
import argparse
import csv
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BINS = [("no orthogroup", None, None), ("0-15", 0, 15), ("16-30", 16, 30), ("31-100", 31, 100),
        ("101-500", 101, 500), (">500", 501, 10 ** 9)]
INK, MUTED, GRID = "#222222", "#6b6b6b", "#e6e6e6"
BAR, LOW = "#2a6f97", "#c0392b"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--explorer-dir", required=True)
    ap.add_argument("--og-table", required=True)
    ap.add_argument("--label", default="Toxo")
    ap.add_argument("--cut", type=int, default=15, help="line drawn at this ortholog count")
    ap.add_argument("--grid-min-strains", type=int, default=10,
                    help="the grid keeps only pangenes with a gene in at least this many strains, with an orthogroup "
                         "(0 = all pangenes, file without the _Nplus_strains suffix)")
    ap.add_argument("--xmax", type=int, default=60, help="panel A shows pangenes with up to this many orthologs")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    cnt = {}
    for r in csv.DictReader(open(a.og_table), delimiter="\t"):
        if r.get("Ortholog Group") and r["Ortholog count"].isdigit():
            cnt[r["Ortholog Group"]] = max(cnt.get(r["Ortholog Group"], 0), int(r["Ortholog count"]))
    fields = json.load(open(os.path.join(a.explorer_dir, "data", "manifest.json")))["gene_fields"]
    qi, pi = fields.index("quality_class"), fields.index("pseudogene")

    rows = []
    for f in sorted(glob.glob(os.path.join(a.explorer_dir, "data", "chunk_*.json"))):
        for p in json.load(open(f)):
            ogs = [o[0] for o in p.get("orthogroups", []) if o[0] in cnt]
            n_orth = max(cnt[o] for o in ogs) if ogs else None
            pg = p.get("panel_genes") or []
            ps = float(pg[0][1]) if pg and pg[0][1] != "" else None
            nv = ni = 0
            pseudo = False
            for recs in p["genes"].values():
                cls = [r[qi] for r in recs]
                pseudo = pseudo or any(r[pi] == "True" for r in recs)
                if "High" in cls:
                    nv += 1
                elif all(c == "Low" for c in cls):
                    ni += 1
            share = ni / (nv + ni) if nv + ni else None
            rows.append(dict(pangene=p["id"], orthogroups=";".join(ogs), orthologs=n_orth, psauron=ps,
                             strains=len(p["genes"]), share_invalid=share,
                             low_confidence="yes" if share is not None and share > 0.5 and not pseudo else "no"))

    def bin_of(r):
        for name, lo, hi in BINS:
            if (lo is None and r["orthologs"] is None) or (lo is not None and r["orthologs"] is not None
                                                           and lo <= r["orthologs"] <= hi):
                return name

    for r in rows:
        r["bin"] = bin_of(r)
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "orthologs_vs_quality.tsv"), "w") as fh:
        keys = ["pangene", "orthogroups", "orthologs", "bin", "psauron", "share_invalid", "low_confidence", "strains"]
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[k] is None else (f"{r[k]:.3f}" if isinstance(r[k], float) else str(r[k]))
                               for k in keys) + "\n")

    names = [b[0] for b in BINS]
    grp = {n: [r for r in rows if r["bin"] == n] for n in names}
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
    for x in ax.flat:
        x.spines[["top", "right"]].set_visible(False)
        x.grid(axis="y", color=GRID, lw=0.8)
        x.set_axisbelow(True)
        x.tick_params(colors=MUTED)

    # A: scatter PSAURON vs orthologs (linear, up to --xmax; jittered so equal counts do not hide each other),
    # coloured by low confidence, with the median PSAURON of each ortholog count
    x0 = ax[0, 0]
    rng = np.random.default_rng(1)
    pts = [r for r in rows if r["orthologs"] is not None and r["psauron"] is not None]
    shown = [r for r in pts if r["orthologs"] <= a.xmax]
    for lc, col, lab, z in (("no", BAR, "other pangenes", 1), ("yes", LOW, "low confidence (invalid ORF in > 50%)", 2)):
        s = [r for r in shown if r["low_confidence"] == lc]
        x0.scatter([r["orthologs"] + rng.uniform(-0.35, 0.35) for r in s], [r["psauron"] for r in s], s=9,
                   alpha=0.35 if lc == "no" else 0.8, color=col, lw=0, label=f"{lab} (n={len(s)})", zorder=z)
    xs = sorted({r["orthologs"] for r in shown})
    med = [np.median([r["psauron"] for r in shown if r["orthologs"] == v]) for v in xs]
    x0.plot(xs, med, color=INK, lw=2, label="median PSAURON per ortholog count", zorder=3)
    x0.axvline(a.cut + 0.5, color=MUTED, ls="--", lw=1)
    x0.text(a.cut + 0.5, 1.04, f" {a.cut} orthologs", color=MUTED, fontsize=8, va="bottom")
    x0.set_xlim(-1, a.xmax + 1)
    x0.text(a.xmax, -0.12, f"{len(pts) - len(shown):,} pangenes with more than {a.xmax} orthologs not shown (see B-D)",
            ha="right", va="top", fontsize=7.5, color=MUTED, transform=x0.get_xaxis_transform())
    x0.set_xlabel("orthologs in VEuPathDB (OrthoMCL)", color=INK)
    x0.set_ylabel("PSAURON of the panel gene", color=INK)
    x0.set_ylim(-0.03, 1.08)
    x0.legend(frameon=False, fontsize=8, loc="lower right")
    x0.set_title("A. PSAURON vs number of orthologs", loc="left", color=INK, fontsize=11)

    # B: PSAURON by ortholog bin (box)
    def boxes(x, key, title, ylab, keep=None):
        data = [[r[key] for r in grp[n] if r[key] is not None and (keep is None or keep(r))] for n in names]
        bp = x.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                       medianprops=dict(color=INK, lw=1.5), whiskerprops=dict(color=MUTED),
                       capprops=dict(color=MUTED), boxprops=dict(edgecolor=MUTED))
        for b in bp["boxes"]:
            b.set_facecolor("#cfe1ec")
        x.set_xticks(range(1, len(names) + 1),
                     [f"{n}\n(n={len(d)}" + (f", {100 * len(d) / len(grp[n]):.0f}%)" if keep and grp[n] else ")")
                      for n, d in zip(names, data)], fontsize=8)
        x.set_xlabel("orthologs in VEuPathDB", color=INK)
        x.set_ylabel(ylab, color=INK)
        x.set_title(title, loc="left", color=INK, fontsize=11)
        for i, d in enumerate(data, 1):
            if d:
                x.text(i, 1.04, f"median {np.median(d):.2f}", ha="center", fontsize=7.5, color=MUTED)
        x.set_ylim(-0.03, 1.1)

    boxes(ax[0, 1], "psauron", "B. PSAURON by number of orthologs", "PSAURON of the panel gene")

    # C: % low PSAURON and % low confidence per bin
    x2 = ax[1, 0]
    pos = np.arange(len(names))
    lowps = [100 * np.mean([r["psauron"] < 0.5 for r in grp[n] if r["psauron"] is not None] or [0]) for n in names]
    lowc = [100 * np.mean([r["low_confidence"] == "yes" for r in grp[n]] or [0]) for n in names]
    w = 0.38
    b1 = x2.bar(pos - w / 2 - 0.01, lowps, w, color=BAR, label="PSAURON < 0.5")
    b2 = x2.bar(pos + w / 2 + 0.01, lowc, w, color=LOW, label="low confidence (invalid ORF in > 50% of strains)")
    for bars in (b1, b2):
        for b in bars:
            x2.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.5, f"{b.get_height():.0f}%", ha="center",
                    fontsize=7.5, color=INK)
    x2.set_xticks(pos, [f"{n}\n(n={len(grp[n])})" for n in names], fontsize=8)
    x2.set_xlabel("orthologs in VEuPathDB", color=INK)
    x2.set_ylabel("% of the pangenes", color=INK)
    x2.legend(frameon=False, fontsize=8, loc="upper right")
    x2.set_title("C. Doubtful pangenes by number of orthologs", loc="left", color=INK, fontsize=11)

    # D: like B, for the share of strains with an invalid ORF; only pangenes with an invalid ORF in 1+ strain
    # (about 90% have none, so with all pangenes every box sits at 0)
    boxes(ax[1, 1], "share_invalid", "D. Share of strains with an invalid ORF (pangenes with 1+ invalid strain)",
          "strains invalid / strains judged", keep=lambda r: r["share_invalid"] is not None and r["share_invalid"] > 0)

    fig.suptitle(f"{a.label}: orthologs (VEuPathDB OrthoMCL, reference gene) vs gene quality, {len(rows):,} pangenes",
                 x=0.01, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(a.out_dir, f"{a.label}_orthologs_vs_quality.png")
    fig.savefig(out, dpi=150)
    print(out)
    grid(rows, names, a)
    scatter3(rows, a)
    for n in names:
        g = grp[n]
        ps = [r["psauron"] for r in g if r["psauron"] is not None]
        print(f"{n}\tpangenes {len(g)}\tmedian_psauron {np.median(ps) if ps else float('nan'):.2f}\t"
              f"psauron<0.5 {sum(1 for v in ps if v < 0.5)}\tlow_conf {sum(1 for r in g if r['low_confidence'] == 'yes')}\t"
              f"median_strains {np.median([r['strains'] for r in g]) if g else 0:.0f}")


INV_BINS = [("none", lambda v: v == 0), ("1-25%", lambda v: 0 < v <= 0.25), ("26-50%", lambda v: 0.25 < v <= 0.5),
            ("51-99%", lambda v: 0.5 < v < 1), ("all (100%)", lambda v: v == 1)]


def grid(rows, names, a):
    """orthologs (rows) x share of strains with an invalid ORF (columns); colour = % of the cell with PSAURON < 0.5.
    Only pangenes with a gene in at least --grid-min-strains strains and an orthogroup."""
    m = a.grid_min_strains
    if m:
        rows = [r for r in rows if r["strains"] >= m and r["orthologs"] is not None]
        names = [n for n in names if n != "no orthogroup"]
    inv = [b[0] for b in INV_BINS]
    cells = {}
    for r in rows:
        if r["share_invalid"] is None:
            continue
        c = next(n for n, t in INV_BINS if t(r["share_invalid"]))
        cells.setdefault((r["bin"], c), []).append(r)
    pct = np.full((len(names), len(inv)), np.nan)
    fig, x = plt.subplots(figsize=(11, 6.8))
    cmap = matplotlib.colormaps["Reds"]
    for i, n in enumerate(names):
        for j, c in enumerate(inv):
            rs = cells.get((n, c), [])
            ps = [r["psauron"] for r in rs if r["psauron"] is not None]
            if ps:
                pct[i, j] = 100 * sum(1 for v in ps if v < 0.5) / len(ps)
    x.imshow(np.nan_to_num(pct, nan=0), cmap=cmap, vmin=0, vmax=100, aspect="auto")
    for i, n in enumerate(names):
        for j, c in enumerate(inv):
            rs = cells.get((n, c), [])
            if not rs:
                x.text(j, i, "-", ha="center", va="center", color=MUTED)
                continue
            ps = [r["psauron"] for r in rs if r["psauron"] is not None]
            dark = not np.isnan(pct[i, j]) and pct[i, j] > 55
            col = "white" if dark else INK
            txt = f"{len(rs):,} pangenes"
            if ps:
                txt += f"\nmedian PSAURON {np.median(ps):.2f}\n{pct[i, j]:.0f}% below 0.5"
            x.text(j, i, txt, ha="center", va="center", fontsize=8.5, color=col, linespacing=1.4)
    x.set_xticks(range(len(inv)), inv)
    x.set_yticks(range(len(names)), [f"{n} (n={sum(1 for r in rows if r['bin'] == n):,})" for n in names])
    x.set_xticks(np.arange(-0.5, len(inv)), minor=True)
    x.set_yticks(np.arange(-0.5, len(names)), minor=True)
    x.grid(which="minor", color="white", lw=2)
    x.tick_params(which="both", length=0, colors=INK)
    for sp in x.spines.values():
        sp.set_visible(False)
    x.xaxis.set_label_position("top")
    x.xaxis.tick_top()
    x.set_xlabel("share of strains with an invalid ORF (most strains invalid = possibly a false gene)", color=INK,
                 labelpad=10)
    x.set_ylabel("orthologs in VEuPathDB (OrthoMCL)", color=INK)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(0, 100))
    cb = fig.colorbar(sm, ax=x, fraction=0.035, pad=0.02)
    cb.set_label("% of the pangenes in the cell with PSAURON < 0.5", color=INK)
    cb.outline.set_visible(False)
    fig.suptitle(f"{a.label}: orthologs x invalid ORFs x PSAURON ({len(rows):,} pangenes"
                 + (f" with a gene in {m}+ strains)" if m else ")"), x=0.01, ha="left",
                 fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(a.out_dir, f"{a.label}_orthologs_invalid_psauron_grid" + (f"_{m}plus_strains" if m else "") + ".png")
    fig.savefig(out, dpi=150)
    print(out)


def spearman_partial(x, y, z):
    """Spearman correlation of x and y once z is accounted for (Pearson of the rank residuals)"""
    from scipy.stats import rankdata
    rx, ry, rz = (rankdata(v) for v in (x, y, z))
    res = lambda r: r - np.polyval(np.polyfit(rz, r, 1), rz)
    return float(np.corrcoef(res(rx), res(ry))[0, 1])


def scatter3(rows, a):
    """the three variables at once: orthologs (x), PSAURON (y), share of invalid strains (colour)"""
    from scipy.stats import spearmanr
    m = a.grid_min_strains
    rs = [r for r in rows if r["orthologs"] is not None and r["psauron"] is not None and r["share_invalid"] is not None
          and r["strains"] >= m]
    o = np.array([r["orthologs"] for r in rs], float)
    p = np.array([r["psauron"] for r in rs])
    v = np.array([r["share_invalid"] for r in rs])
    rng = np.random.default_rng(1)
    cap = a.xmax
    xpos = np.where(o <= cap, o + rng.uniform(-0.35, 0.35, len(o)), cap + 3 + rng.uniform(-2, 2, len(o)))

    fig, x = plt.subplots(figsize=(12, 6.8))
    x.spines[["top", "right"]].set_visible(False)
    x.grid(axis="y", color=GRID, lw=0.8)
    x.set_axisbelow(True)
    zero = v == 0
    x.scatter(xpos[zero], p[zero], s=8, color="#b9c4cc", alpha=0.45, lw=0, label=f"no invalid strain (n={zero.sum():,})")
    order = np.argsort(v[~zero])
    sc = x.scatter(xpos[~zero][order], p[~zero][order], s=16, c=v[~zero][order], cmap="Reds", vmin=0, vmax=1,
                   edgecolors="#555555", linewidths=0.3, label=f"1+ invalid strain (n={(~zero).sum():,})")
    cb = fig.colorbar(sc, ax=x, fraction=0.035, pad=0.02)
    cb.set_label("share of strains with an invalid ORF", color=INK)
    cb.outline.set_visible(False)
    x.axvline(a.cut + 0.5, color=MUTED, ls="--", lw=1)
    x.text(a.cut + 0.5, 1.04, f" {a.cut} orthologs", color=MUTED, fontsize=8, va="bottom")
    x.axvline(cap + 0.75, color=GRID, lw=6)
    ticks = list(range(0, cap + 1, 10))
    x.set_xticks(ticks + [cap + 3], [str(t) for t in ticks] + [f">{cap}\n(n={(o > cap).sum():,})"])
    x.set_xlim(-1, cap + 6)
    x.set_ylim(-0.03, 1.1)
    x.set_xlabel("orthologs in VEuPathDB (OrthoMCL)", color=INK)
    x.set_ylabel("PSAURON of the panel gene", color=INK)
    x.legend(frameon=False, fontsize=8, loc="lower right")

    r_op = spearmanr(o, p)[0]
    r_ov = spearmanr(o, v)[0]
    r_pv = spearmanr(p, v)[0]
    txt = ("Spearman correlation (all points, ranks)\n"
           f"  orthologs x PSAURON        {r_op:+.2f}   (without invalid share: {spearman_partial(o, p, v):+.2f})\n"
           f"  orthologs x invalid share  {r_ov:+.2f}   (without PSAURON: {spearman_partial(o, v, p):+.2f})\n"
           f"  PSAURON x invalid share    {r_pv:+.2f}   (without orthologs: {spearman_partial(p, v, o):+.2f})")
    x.text(0.99, 0.30, txt, transform=x.transAxes, ha="right", va="bottom", multialignment="left", fontsize=8, family="monospace",
           color=INK, bbox=dict(facecolor="white", edgecolor=GRID, boxstyle="round,pad=0.5"))
    fig.suptitle(f"{a.label}: orthologs, PSAURON and invalid ORFs, one point per pangene ({len(rs):,} pangenes"
                 + (f" with a gene in {m}+ strains)" if m else ")"), x=0.01, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(a.out_dir, f"{a.label}_orthologs_psauron_invalid_scatter" + (f"_{m}plus_strains" if m else "")
                       + ".png")
    fig.savefig(out, dpi=150)
    print(out)
    print(txt)


if __name__ == "__main__":
    main()
