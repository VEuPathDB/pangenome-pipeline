#!/usr/bin/env python3
"""Produce a 6-panel comparative QC figure from assess_assemblies.py output.

Subpanels:
  (a) Genome size & contig count (dual axis)
  (b) Contig N50
  (c) GC content (%)
  (d) Gene / mRNA / CDS counts
  (e) Exons & introns (derived)
  (f) Protein length distribution (boxplot)
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Generate 6-panel comparative assembly QC figure."
    )
    parser.add_argument(
        "--metrics",
        default="assembly_metrics.tsv",
        help="Path to assembly_metrics.tsv from assess_assemblies.py (default: assembly_metrics.tsv)",
    )
    parser.add_argument(
        "--proteins",
        default="protein_lengths.tsv",
        help="Path to protein_lengths.tsv from assess_assemblies.py (default: protein_lengths.tsv)",
    )
    parser.add_argument(
        "--output",
        default="assembly_assessment.png",
        help="Output image path (default: assembly_assessment.png)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=200,
        help="Image resolution DPI (default: 200)",
    )
    args = parser.parse_args()

    if not os.path.exists(args.metrics):
        print(f"Error: Metrics file not found: {args.metrics}", file=sys.stderr)
        sys.exit(1)

    m = pd.read_csv(args.metrics, sep="\t")
    if m.empty:
        print(f"Error: Metrics table {args.metrics} is empty.", file=sys.stderr)
        sys.exit(1)

    has_proteins = os.path.exists(args.proteins) and os.path.getsize(args.proteins) > 0
    if has_proteins:
        pl = pd.read_csv(args.proteins, sep="\t")
    else:
        pl = pd.DataFrame(columns=["accession", "protein_length"])

    labels = m["accession"].astype(str).tolist()
    n_genomes = len(labels)
    x = np.arange(n_genomes)

    # Dynamic figure width based on number of assemblies
    fig_width = max(16, min(36, 10 + n_genomes * 1.2))
    fig, axes = plt.subplots(2, 3, figsize=(fig_width, 11))

    # --- (a) Genome size + contig count (dual axis) ---
    ax = axes[0, 0]
    bar_w = 0.35
    b1 = ax.bar(
        x - bar_w / 2,
        m["genome_size_bp"] / 1e6,
        bar_w,
        label="Size (Mb)",
        color="steelblue",
        edgecolor="none",
    )
    ax2 = ax.twinx()
    b2 = ax2.bar(
        x + bar_w / 2,
        m["num_contigs"],
        bar_w,
        label="Contigs",
        color="orange",
        alpha=0.85,
        edgecolor="none",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Genome size (Mb)", color="steelblue", fontweight="bold")
    ax2.set_ylabel("Contig count", color="darkorange", fontweight="bold")
    ax.set_title("A. Genome size & contig count", fontweight="bold")
    ax.legend(handles=[b1, b2], loc="upper left", frameon=True)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # --- (b) Contig N50 ---
    ax = axes[0, 1]
    ax.bar(x, m["n50"] / 1e6, color="seagreen", width=0.55, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("N50 (Mb)", fontweight="bold")
    ax.set_title("B. Contig N50", fontweight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # --- (c) GC content ---
    ax = axes[0, 2]
    ax.bar(x, m["gc_pct"], color="mediumpurple", width=0.55, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("GC content (%)", fontweight="bold")
    ax.set_ylim(0, 100)
    ax.set_title("C. GC content", fontweight="bold")
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # --- (d) Feature counts (gene, mRNA, CDS) ---
    ax = axes[1, 0]
    w = 0.25
    ax.bar(x - w, m["num_genes"], w, label="Genes", color="#4C72B0", edgecolor="none")
    ax.bar(x, m["num_mRNAs"], w, label="mRNAs", color="#55A868", edgecolor="none")
    ax.bar(x + w, m["num_CDS"], w, label="CDS", color="#C44E52", edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Feature count", fontweight="bold")
    ax.set_title("D. Gene / mRNA / CDS counts", fontweight="bold")
    ax.legend(frameon=True)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # --- (e) Exons & Introns ---
    ax = axes[1, 1]
    w = 0.35
    ax.bar(x - w / 2, m["num_exons"], w, label="Exons", color="teal", edgecolor="none")
    ax.bar(
        x + w / 2,
        m["num_introns_derived"],
        w,
        label="Introns (derived)",
        color="coral",
        edgecolor="none",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_ylabel("Count", fontweight="bold")
    ax.set_title("E. Exons & derived introns", fontweight="bold")
    ax.legend(frameon=True)
    ax.grid(axis="y", linestyle=":", alpha=0.4)

    # --- (f) Protein length distribution ---
    ax = axes[1, 2]
    if has_proteins and not pl.empty:
        data = [
            pl.loc[pl["accession"].astype(str) == a, "protein_length"].dropna().values
            for a in labels
        ]
        # Only plot if at least one accession has lengths
        if any(len(d) > 0 for d in data):
            ax.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True)
            ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
            ax.set_ylabel("Protein length (aa)", fontweight="bold")
            ax.set_title("F. Protein length distribution", fontweight="bold")
            ax.grid(axis="y", linestyle=":", alpha=0.4)
        else:
            ax.text(0.5, 0.5, "No protein sequences found", ha="center", va="center")
            ax.set_title("F. Protein length distribution", fontweight="bold")
    else:
        ax.text(0.5, 0.5, "Protein lengths file not provided", ha="center", va="center")
        ax.set_title("F. Protein length distribution", fontweight="bold")

    plt.tight_layout()
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    plt.savefig(args.output, dpi=args.dpi)
    print(f"Successfully generated figure: {args.output}")


if __name__ == "__main__":
    main()
