#!/usr/bin/env python3
"""Gene properties of the pangenes by pangene occupancy (domain evidence and protein length).

Only pangenes with exactly one gene in every genome they are present in are used, so each pangene of
occupancy k has exactly k genes (pangenes with more than one gene in a genome are left out).

  B  Proportion of genes with an InterPro/PFAM domain - one bar per occupancy (mean over the
     pangenes of that occupancy of "genes with a domain / genes in the pangene").
  C  Protein length in amino acids - one boxplot per occupancy. Each pangene contributes ONE value:
     the median of its genes' protein lengths (a gene's protein length = summed CDS length / 3, of its
     longest transcript).

Bars and boxes are coloured by category (Core / Soft-core / Shell / Cloud, same definition as
plot_occupancy.py). The Tukey HSD p-value of Core versus Cloud is written on each panel; note that the
genes of one pangene come from different genomes and are not independent samples.

Inputs:
  --genes-matrix  pangene_matrix_genes.tr.tab from get_pangenes (cells = gene IDs, "-" = absent)
  --gff-dir       GFFs whose gene IDs match the matrix, named <genome><--gff-suffix>
Writes <out-prefix>.png/.svg and <out-prefix>.tsv (one row per pangene used).

Needs numpy, scipy, matplotlib.
"""

import argparse
import csv
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter, NullFormatter
from matplotlib.transforms import blended_transform_factory
from plot_occupancy import COLORS, ORDER, categorize
from scipy.stats import tukey_hsd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.gff import DBXREF_RE, ID_RE, PARENT_RE


def parse_gff(path):
    """gene ID -> (has_domain, protein length in aa) from the CDS lines (gene -> mRNA -> CDS)."""
    mrna_gene, cds_bp, domain = {}, {}, set()
    genes = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = f[8]
            if f[2] == "gene":
                i = ID_RE.search(attrs)
                if i:
                    genes.add(i.group(1))
            elif f[2] in ("mRNA", "transcript"):
                i, p = ID_RE.search(attrs), PARENT_RE.search(attrs)
                if i and p:
                    mrna_gene[i.group(1)] = p.group(1)
            elif f[2] == "CDS":
                p = PARENT_RE.search(attrs)
                if not p:
                    continue
                m = p.group(1)
                cds_bp[m] = cds_bp.get(m, 0) + int(f[4]) - int(f[3]) + 1
                d = DBXREF_RE.search(attrs)
                if d and any(t.startswith(("InterPro:", "PFAM:")) for t in d.group(1).split(",")):
                    domain.add(m)
    out = {g: [False, 0] for g in genes}
    for m, bp in cds_bp.items():
        g = mrna_gene.get(m)
        if g in out:
            out[g][1] = max(out[g][1], bp / 3.0)
            if m in domain:
                out[g][0] = True
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument(
        "--gff-suffix",
        default="_for_get_pangenes.gff",
        help="GFF file name after the genome name (default _for_get_pangenes.gff, the names in the get_pangenes input folder)",
    )
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument(
        "--core-missing", type=int, default=0, help="core = missing from at most this many genomes (default 0)"
    )
    ap.add_argument(
        "--softcore-missing", type=int, default=1, help="soft-core = missing from at most this many genomes (default 1)"
    )
    ap.add_argument("--cloud-max", type=int, default=2, help="cloud = present in at most this many genomes (default 2)")
    ap.add_argument(
        "--yscale",
        choices=["linear", "log"],
        default="linear",
        help="y-axis of the protein-length panel (default linear; on linear the axis stops at the 99th percentile and the number of pangenes above it is written on the panel)",
    )
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    with open(args.genes_matrix, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    cols = [(i, g.replace("_for_get_pangenes", "")) for i, g in enumerate(rows[0]) if i > 0 and g.strip()]
    n_genomes = len(cols)

    gene_data = {}
    for _, genome in cols:
        path = os.path.join(args.gff_dir, genome + args.gff_suffix)
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")
        gene_data[genome] = parse_gff(path)
        print(f"{genome}: {len(gene_data[genome])} genes", file=sys.stderr)

    used, multi, missing = [], 0, 0
    for r in rows[1:]:
        if r[0].startswith("chr"):
            continue
        members = []
        for i, genome in cols:
            cell = r[i] if i < len(r) else ""
            if cell in ("", "-", "NA"):
                continue
            ids = cell.split(",")
            members.append((genome, ids))
        if not members:
            continue
        if any(len(ids) > 1 for _, ids in members):
            multi += 1
            continue
        vals = [gene_data[g].get(ids[0]) for g, ids in members]
        if any(v is None for v in vals):
            missing += 1
            continue
        occ = len(members)
        dom = float(np.mean([v[0] for v in vals]))
        length = float(np.median([v[1] for v in vals]))
        used.append(
            (
                r[0],
                occ,
                categorize(occ, n_genomes, args.core_missing, args.softcore_missing, args.cloud_max),
                dom,
                length,
            )
        )

    print(
        f"{len(used)} pangenes used; {multi} left out (more than one gene in a genome)"
        + (f"; {missing} skipped (gene not found in the GFFs)" if missing else ""),
        file=sys.stderr,
    )
    if not used:
        sys.exit("ERROR: no pangene left; check that --gff-dir IDs match the matrix gene IDs")

    with open(args.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pangene", "occupancy", "category", "domain_proportion", "median_protein_aa"])
        for u in used:
            w.writerow([u[0], u[1], u[2], f"{u[3]:.3f}", f"{u[4]:.1f}"])

    occ_all = np.array([u[1] for u in used])
    cat_all = np.array([u[2] for u in used])
    dom_all = np.array([u[3] for u in used])
    len_all = np.array([u[4] for u in used])
    ks = [k for k in range(1, n_genomes + 1) if (occ_all == k).any()]
    cat_of = {k: cat_all[occ_all == k][0] for k in ks}
    counts = [int((occ_all == k).sum()) for k in ks]

    def core_vs_cloud(values):
        groups = [values[cat_all == c] for c in ORDER if (cat_all == c).sum() > 1]
        names = [c for c in ORDER if (cat_all == c).sum() > 1]
        if "Core" not in names or "Cloud" not in names:
            return None
        res = tukey_hsd(*groups)
        return float(res.pvalue[names.index("Core"), names.index("Cloud")])

    def fmt_p(p):
        return "n/a" if p is None else ("p < 0.001" if p < 0.001 else f"p = {p:.3f}")

    fig, (axb, axc) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    xs = np.arange(len(ks))
    cloud_x = [x for x, k in zip(xs, ks) if cat_of[k] == "Cloud"]
    core_x = [x for x, k in zip(xs, ks) if cat_of[k] == "Core"]

    def significance_line(ax, p):
        """Bracket from the cloud classes to the core class, with the Tukey HSD p-value above it."""
        if p is None or not cloud_x or not core_x:
            return
        x0, x1 = float(np.mean(cloud_x)), float(np.mean(core_x))
        ax.plot(
            [x0, x0, x1, x1],
            [0.90, 0.93, 0.93, 0.90],
            transform=blended_transform_factory(ax.transData, ax.transAxes),
            color="black",
            lw=1,
            clip_on=False,
        )
        ax.text(
            (x0 + x1) / 2,
            0.94,
            fmt_p(p),
            transform=blended_transform_factory(ax.transData, ax.transAxes),
            ha="center",
            va="bottom",
            fontsize=9,
        )

    axb.bar(xs, [dom_all[occ_all == k].mean() for k in ks], width=0.8, color=[COLORS[cat_of[k]] for k in ks])
    axb.set_ylabel("Proportion of genes with\na domain", fontsize=12)
    axb.set_ylim(0, 1.2)
    axb.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    significance_line(axb, core_vs_cloud(dom_all))

    bp = axc.boxplot(
        [len_all[occ_all == k] for k in ks],
        positions=xs,
        widths=0.7,
        patch_artist=True,
        medianprops={"color": "black", "lw": 1.5},
        flierprops={"marker": "o", "ms": 2.5, "alpha": 0.4, "mec": "none", "mfc": "0.3"},
    )
    for patch, k in zip(bp["boxes"], ks):
        patch.set_facecolor(COLORS[cat_of[k]])
        patch.set_alpha(0.9)
    plain = FuncFormatter(lambda v, _: f"{v:,.0f}")
    if args.yscale == "log":
        axc.set_yscale("log")
        axc.yaxis.set_major_formatter(plain)
        axc.yaxis.set_minor_formatter(NullFormatter())
        axc.set_ylabel("Protein length (aa, log scale),\nmedian per pangene", fontsize=12)
        axc.set_ylim(top=len_all.max() * 3)
    else:
        top = float(np.percentile(len_all, 99))
        axc.set_ylim(0, top / 0.85)
        axc.yaxis.set_major_formatter(plain)
        for fl in bp["fliers"]:
            xd, yd = fl.get_xdata(), fl.get_ydata()
            keep = yd <= top
            fl.set_data(xd[keep], yd[keep])
        n_out = int((len_all > top).sum())
        axc.set_ylabel("Protein length (aa),\nmedian per pangene", fontsize=12)
        axc.text(
            0.5,
            0.80,
            f"{n_out} pangenes above {top:,.0f} aa not shown",
            transform=axc.transAxes,
            ha="center",
            fontsize=8,
        )
    significance_line(axc, core_vs_cloud(len_all))

    axc.set_xticks(xs)
    axc.set_xticklabels([f"{k}\n(n={c})" for k, c in zip(ks, counts)], fontsize=8)
    axc.set_xlabel("Pangene occupancy", fontsize=12)
    for ax in (axb, axc):
        ax.spines[["top", "right"]].set_visible(False)
    present = [c for c in ORDER if c in cat_of.values()]
    axb.legend(
        [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in present],
        present,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=len(present),
        frameon=False,
    )
    if args.title:
        fig.suptitle(args.title, x=0.01, ha="left", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
