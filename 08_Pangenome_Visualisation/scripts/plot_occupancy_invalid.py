#!/usr/bin/env python3
"""Pangenes with mostly invalid genes, by pangene occupancy.

A gene is "invalid" when it is pseudogene-like, by either mark (same definition as
07_Post_GET_PANGENES_adjusting/scripts/explore_pseudogenes.py):
  biotype     ebi_biotype=pseudogene or gene_biotype=pseudogene
  broken ORF  valid_ORFs=0 (or valid_ORF=False), a Liftoff transfer with no valid open reading frame
A gene can only be judged if it has one of the two marks or an ORF result. Genes taken straight from a
genome's own annotation (no Liftoff) have no ORF result, so they count as "not judged" and stay out of the
denominator, unless their biotype says pseudogene.

For each pangene: share of invalid genes among the judged genes of all its genomes (every gene counts,
several genes in one genome included). A pangene is "mostly invalid" when that share is above --min-share
(default 0.5). Pangenes with no judged gene are left out and reported.

--mark overlap uses another definition of a bad gene, for the tagged (7b) GFFs: a gene is "unsupported" when
it is NOT from the reference layer (source `main`), has no evidence and overlaps another gene (`intersect=`
is not NA). Evidence = an InterPro/Pfam entry in the CDS Dbxref OR a specific description (anything other than
empty, "hypothetical protein" and the other generic texts of filter_final_annotation.py). Main genes, and
genes with evidence, count as good even when they overlap. Every gene can be judged this way.

One bar per occupancy (1..N genomes), coloured by category (Core / Soft-core / Shell / Cloud, same
definition as plot_occupancy.py); the hatched part of each bar is the pangenes that are mostly invalid.

Inputs:
  --genes-matrix  pangene_matrix_genes.tr.tab from get_pangenes (cells = gene IDs, "-" = absent)
  --gff-dir       GFFs whose gene IDs match the matrix, named <genome><--gff-suffix>
--single-copy leaves out the pangenes that have more than one gene in some genome.
Writes <out-prefix>.png/.svg and <out-prefix>.tsv (occupancy, category, pangenes, mostly_invalid) and
<out-prefix>_pangenes.tsv (one row per pangene).

Needs numpy, matplotlib.
"""
import argparse
import csv
import os
import re
import sys
from urllib.parse import unquote_plus

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_occupancy import categorize, COLORS, ORDER

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")


def gene_status(path):
    """gene ID -> True (invalid), False (valid) or None (not judged)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = f[8]
            i = ID_RE.search(a)
            if not i:
                continue
            biotype = "ebi_biotype=pseudogene" in a or "gene_biotype=pseudogene" in a
            known = broken = False
            m = re.search(r"valid_ORFs=(\d+)", a)
            if m:
                known, broken = True, int(m.group(1)) == 0
            if "valid_ORF=False" in a:
                known, broken = True, True
            elif "valid_ORF=True" in a:
                known = True
            out[i.group(1)] = True if (biotype or broken) else (False if known else None)
    return out


PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
INTERSECT_RE = re.compile(r"(?:^|;)intersect=([^;]+)")
DESCRIPTION_RE = re.compile(r"(?:^|;)description=([^;]+)")
PRODUCT_RE = re.compile(r"(?:^|;)product=([^;]+)")
# Same list as filter_final_annotation.py: a description this generic counts as no evidence.
GENERIC_DESCRIPTIONS = {
    "hypothetical protein", "hypothetical protein, conserved", "conserved hypothetical protein",
    "unspecified product", "unknown function", "unnamed protein product", "",
}


def gene_status_overlap(path):
    """gene ID -> True (not main, no domain, no specific description, overlaps another gene) or False; untagged genes are left out."""
    genes, mrna_gene, with_domain = {}, {}, set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            i = ID_RE.search(f[8])
            if f[2] == "gene" and i:
                m = INTERSECT_RE.search(f[8])
                if m:   # genes without the tag (e.g. non-coding RNA genes typed gene) are not judged
                    d = DESCRIPTION_RE.search(f[8]) or PRODUCT_RE.search(f[8])
                    specific = (unquote_plus(d.group(1).strip().lower()) if d else "") not in GENERIC_DESCRIPTIONS
                    genes[i.group(1)] = (f[1] != "main", m.group(1) != "NA", specific)
            elif f[2] in ("mRNA", "transcript") and i:
                p = PARENT_RE.search(f[8])
                if p:
                    mrna_gene[i.group(1)] = p.group(1)
            elif f[2] == "CDS":
                p = PARENT_RE.search(f[8])
                if p and ("InterPro:" in f[8] or "PFAM:" in f[8] or "Pfam:" in f[8]):
                    with_domain.add(p.group(1))
    if not genes:
        sys.exit(f"ERROR: no gene with an intersect= tag in {path}; --mark overlap needs the tagged (7b) GFFs")
    domain_genes = {mrna_gene.get(m, m) for m in with_domain}
    return {g: (not_main and overlaps and not specific and g not in domain_genes)
            for g, (not_main, overlaps, specific) in genes.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_for_get_pangenes.gff")
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument("--min-share", type=float, default=0.5,
                    help="a pangene is mostly invalid when MORE than this share of its judged genes are invalid (default 0.5)")
    ap.add_argument("--core-missing", type=int, default=0, help="core = missing from at most this many genomes (default 0)")
    ap.add_argument("--softcore-missing", type=int, default=1, help="soft-core = missing from at most this many genomes (default 1)")
    ap.add_argument("--cloud-max", type=int, default=2, help="cloud = present in at most this many genomes (default 2)")
    ap.add_argument("--style", choices=["counts", "percent", "invalid", "share"], default="counts",
                    help="counts = bar height is the number of pangenes; percent = every bar goes 0-100 and shows the share that is mostly invalid; invalid = only the mostly invalid pangenes, counted per occupancy; share = one stacked bar per category, split by the share of invalid genes of each pangene")
    ap.add_argument("--single-copy", action="store_true",
                    help="only pangenes with exactly one gene in every genome they are present in (leaves out the multi-gene pangenes)")
    ap.add_argument("--mark", choices=["orf", "overlap"], default="orf",
                    help="orf = pseudogene-like genes (broken ORF or pseudogene biotype); overlap = not main, no domain and overlapping another gene (tagged 7b GFFs)")
    ap.add_argument("--title", default="")
    args = ap.parse_args()
    if args.mark == "overlap":
        what, bad_gene, good_gene = "unsupported genes", "unsupported", "supported"
    else:
        what, bad_gene, good_gene = "invalid ORFs", "invalid", "valid"

    with open(args.genes_matrix, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    cols = [(i, g.replace("_for_get_pangenes", "")) for i, g in enumerate(rows[0]) if i > 0 and g.strip()]
    n = len(cols)
    status = {}
    for _, genome in cols:
        path = os.path.join(args.gff_dir, genome + args.gff_suffix)
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")
        status[genome] = gene_status_overlap(path) if args.mark == "overlap" else gene_status(path)

    total = np.zeros(n + 1, dtype=int)
    bad = np.zeros(n + 1, dtype=int)
    skipped = n_multi = 0
    per_pangene = []
    for r in rows[1:]:
        if r[0].startswith("chr"):
            continue
        occ, genes, multi_gene = 0, [], False
        for i, genome in cols:
            cell = r[i] if i < len(r) else ""
            if cell not in ("", "-", "NA"):
                occ += 1
                ids = cell.split(",")
                multi_gene |= len(ids) > 1
                genes += [status[genome].get(g) for g in ids]
        if not occ:
            continue
        if args.single_copy and multi_gene:
            n_multi += 1
            continue
        judged = [s for s in genes if s is not None]
        if not judged:
            skipped += 1
            continue
        share = sum(judged) / len(judged)
        mostly = share > args.min_share
        total[occ] += 1
        bad[occ] += mostly
        per_pangene.append((r[0], occ, len(genes), len(judged), sum(judged), share, mostly))
    print(f"{total.sum()} pangenes judged; {skipped} left out (no gene could be judged)"
          + (f"; {n_multi} left out (more than one gene in a genome)" if args.single_copy else ""), file=sys.stderr)

    ks = np.arange(1, n + 1)
    cats = [categorize(k, n, args.core_missing, args.softcore_missing, args.cloud_max) for k in ks]
    with open(args.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["occupancy", "category", "pangenes", "mostly_invalid"])
        for k, c in zip(ks, cats):
            w.writerow([k, c, int(total[k]), int(bad[k])])
    with open(args.out_prefix + "_pangenes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["pangene", "occupancy", "genes", "genes_judged", "genes_invalid", "share_invalid", "mostly_invalid"])
        for p in per_pangene:
            w.writerow([p[0], p[1], p[2], p[3], p[4], f"{p[5]:.2f}", int(p[6])])

    print(f"\n  {'category':10s}{'pangenes':>10s}{'mostly invalid':>16s}{'%':>8s}")
    for cat in ORDER:
        sel = [k for k, c in zip(ks, cats) if c == cat]
        t, b = int(sum(total[k] for k in sel)), int(sum(bad[k] for k in sel))
        if t:
            print(f"  {cat:10s}{t:10d}{b:16d}{100 * b / t:7.1f}%")

    if args.style == "share":
        cat_of = {k: c for k, c in zip(ks, cats)}
        pc = np.array([cat_of[p[1]] for p in per_pangene])
        sh = np.array([p[5] for p in per_pangene])
        classes = [(lambda x: x == 0, f"0% (all genes {good_gene})"),
                   (lambda x: (x > 0) & (x <= 0.5), "up to 50%"),
                   (lambda x: (x > 0.5) & (x < 1), "more than 50%, not all"),
                   (lambda x: x == 1, f"100% (all genes {bad_gene})")]
        shades = ["#e8e8e8", "#fdd0a2", "#f16913", "#7f2704"]
        show = [c for c in reversed(ORDER) if (pc == c).any()]   # Cloud left, Core right
        fig, ax = plt.subplots(figsize=(7.5, 5.2))
        bottom = np.zeros(len(show))
        for (fn, name), col in zip(classes, shades):
            share = np.array([100 * fn(sh[pc == c]).mean() for c in show])
            ax.bar(range(len(show)), share, bottom=bottom, color=col, width=0.7, label=name, edgecolor="white", lw=0.5)
            for x, (v, b) in enumerate(zip(share, bottom)):
                if v >= 3:
                    ax.text(x, b + v / 2, f"{v:.0f}%", ha="center", va="center", fontsize=8,
                            color="white" if col in shades[2:] else "black")
            bottom += share
        ax.set_xticks(range(len(show)))
        ax.set_xticklabels([f"{c}\n(n={int((pc == c).sum())})" for c in show])
        ax.set_ylabel("Percentage of pangenes", fontsize=12)
        ax.set_ylim(0, 100)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(title=f"{bad_gene.capitalize()} genes in the pangene", loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
        fig.tight_layout()
        for ext in ("png", "svg"):
            fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
        print(f"\nWrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)
        return

    fig, ax = plt.subplots(figsize=(9, 5))
    colors = [COLORS[c] for c in cats]
    if args.style == "invalid":
        ax.bar(ks, bad[1:], width=0.8, color=colors)
        for k in ks:
            ax.text(k, bad[k], f"{bad[k]}", ha="center", va="bottom", fontsize=8)
        ax.set_ylim(0, bad.max() * 1.08)
        ax.set_xticks(ks)
        ax.set_ylabel(f"Number of pangenes mostly {bad_gene}" if args.mark == "overlap" else f"Number of pangenes with mostly {what}", fontsize=12)
    elif args.style == "percent":
        share = 100 * bad[1:] / np.maximum(total[1:], 1)
        ax.bar(ks, share, width=0.8, color=colors)
        for k, sh in zip(ks, share):
            if total[k]:
                ax.text(k, sh, f"{sh:.0f}%", ha="center", va="bottom", fontsize=8)
        ax.set_ylim(0, 100)
        ax.set_xticks(ks)
        ax.set_xticks(ks)
        ax.set_ylabel(f"Pangenes mostly {bad_gene} (%)" if args.mark == "overlap" else f"Pangenes with mostly {what} (%)", fontsize=12)
    else:
        ax.bar(ks, total[1:] - bad[1:], width=0.8, color=colors)
        ax.bar(ks, bad[1:], bottom=total[1:] - bad[1:], width=0.8, color=colors, hatch="////", edgecolor="white", lw=0)
        for k in ks:
            if bad[k]:
                ax.text(k, total[k], f"{bad[k]}", ha="center", va="bottom", fontsize=8)
        ax.set_ylim(0, total.max() * 1.08)
        ax.set_xticks(ks)
        ax.set_ylabel("Number of pangenes", fontsize=12)
    ax.set_xlabel("Pangene occupancy", fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    present = [c for c in ORDER if c in cats]
    handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[c]) for c in present]
    labels = list(present)
    if args.style == "counts":
        handles.append(plt.Rectangle((0, 0), 1, 1, facecolor="0.6", hatch="////", edgecolor="white", lw=0))
        labels.append(f"> {args.min_share * 100:.0f}% {bad_gene} genes")
    ax.legend(handles, labels,
              loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    if args.title:
        ax.set_title(args.title, loc="left", fontsize=14)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"\nWrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
