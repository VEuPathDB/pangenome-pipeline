#!/usr/bin/env python3
"""
Plots comparing this pipeline's output before and after filter_final_annotation.py (7a) -
run this after 7a, as a sanity check that the filtering did what it should across every
strain, not just the ones checked by hand.

Reuses filter_final_annotation.py's own low-level GFF-attribute helpers directly (imported, not
copy-pasted) - it needs to live in the same folder as filter_final_annotation.py to import them.
Gene-level parsing itself is NOT reused: the "before" file (<genome>_final_annotation.gff) still
types genes `protein_coding_gene`, but 7a's own output is already relabeled to `gene` (so
get_pangenes can cluster it directly) - this script parses each side with its own feature type.

Needs matplotlib (`pip install matplotlib`) - everything else is standard library.

Produces, in --out-dir:
    01_gene_counts.png       total genes, before vs. after, per strain
    02_extra_copies.png      extra-copy genes remaining, before vs. after, per strain
    03_removal_reasons.png   why genes were removed, per strain (stacked)
    04_overlap_severity.png  CDS-overlap severity, before vs. after, per strain
    05_source_composition.png  where the genes came from, before vs. after, per strain
    summary.tsv              the numbers behind every plot above, one strain per row

Usage:
    python3 plot_filter_qc.py \\
        --strains genome_order.txt \\
        --original-dir expanded_annotation \\
        --filtered-dir filtered_annotation \\
        --out-dir filter_qc_plots
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("ERROR: this script needs matplotlib - install it with 'pip install matplotlib'")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from annotation_common import ID_RE, PARENT_RE, find_cds_overlaps, float_attr, int_attr, source_counts, valid_orf

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.io import read_strains

REASON_CATEGORIES = [
    (r"^extra copy with no valid ORF", "extra copy: no valid ORF"),
    (r"^rescued transfer, Liftoff coverage", "rescued: below min coverage"),
    (r"overlaps the reference \(", "overlap: vs. main, no evidence"),
    (r"has evidence;", "overlap: other side has evidence"),
    (r"neither has evidence, \S+ more supported", "overlap: good-transfer/size tiebreak (neither has evidence)"),
]


def categorize_reason(reason):
    for pattern, label in REASON_CATEGORIES:
        if re.search(pattern, reason):
            return label
    return "other"


def parse_gff_for_qc(path, gene_feature_type):
    """Same fields filter_final_annotation.py's own parse_gff produces, just parameterized on
    which feature type marks a gene-level line - "protein_coding_gene" for the raw pre-filter
    file, "gene" for 7a's already-relabeled output."""
    genes = {}
    mrna_to_gene = {}
    cds_by_mrna = defaultdict(list)
    mrna_has_domain = defaultdict(bool)

    with open(path) as fh:
        raw_lines = fh.readlines()

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqname, source, feature, start, end, _score, strand = fields[:7]
        attrs = fields[8]
        m = ID_RE.search(attrs)
        feat_id = m.group(1) if m else None

        if feature == gene_feature_type and feat_id:
            genes[feat_id] = {
                "id": feat_id,
                "seqname": seqname,
                "source": source,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "extra_copy_number": int_attr(attrs, "extra_copy_number"),
                "valid_ORF": valid_orf(attrs),
                "coverage": float_attr(attrs, "coverage"),
            }
        elif feature == "mRNA" and feat_id:
            parent = PARENT_RE.search(attrs)
            if parent:
                mrna_to_gene[feat_id] = parent.group(1)

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        feature, start, end, attrs = fields[2], fields[3], fields[4], fields[8]
        if feature != "CDS":
            continue
        parent = PARENT_RE.search(attrs)
        if not parent:
            continue
        mrna_id = parent.group(1)
        if mrna_to_gene.get(mrna_id) not in genes:
            continue
        cds_by_mrna[mrna_id].append((int(start), int(end)))
        if "InterPro:" in attrs or "PFAM:" in attrs or "Pfam:" in attrs:
            mrna_has_domain[mrna_id] = True

    for gene_id, gene in genes.items():
        gene_mrnas = [m for m, g in mrna_to_gene.items() if g == gene_id]
        best_mrna, best_bp = None, -1
        for mrna_id in gene_mrnas:
            segs = cds_by_mrna.get(mrna_id, [])
            bp = sum(e - s + 1 for s, e in segs)
            if bp > best_bp:
                best_bp, best_mrna = bp, mrna_id
        gene["cds_segments"] = sorted(cds_by_mrna.get(best_mrna, []))
        gene["cds_bp"] = max(0, best_bp)
        gene["has_domain"] = any(mrna_has_domain.get(m, False) for m in gene_mrnas)

    return genes


def severity_counts(genes, min_overlap_bp):
    """Same thresholds as filter_final_annotation.py's overlap-resolution step - minor/
    substantial/contained, sized against the smaller gene's own CDS."""
    pairs = find_cds_overlaps(genes, min_overlap_bp)
    counts = {"minor": 0, "substantial": 0, "contained": 0}
    for gene1_id, gene2_id, ov_bp in pairs:
        g1, g2 = genes[gene1_id], genes[gene2_id]
        smaller_bp = min(g1["cds_bp"], g2["cds_bp"])
        if smaller_bp == 0:
            continue
        pct = ov_bp / smaller_bp
        if pct >= 0.8:
            counts["contained"] += 1
        elif pct >= 0.2:
            counts["substantial"] += 1
        else:
            counts["minor"] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--strains",
        required=True,
        help="Comma-separated genome names or a file with one per line - first one must be the reference",
    )
    parser.add_argument(
        "--original-dir",
        required=True,
        help="Directory with <strain>_final_annotation.gff (this pipeline's raw output, before filtering)",
    )
    parser.add_argument(
        "--filtered-dir",
        required=True,
        help="Directory with <strain>_final_annotation_filtered.gff and <strain>_filter_report.tsv (7a's output)",
    )
    parser.add_argument(
        "--min-overlap-bp",
        type=int,
        default=20,
        help="Same --min-overlap-bp used for the real filter run (default: 20)",
    )
    parser.add_argument("--out-dir", default="./filter_qc_plots", help="Where to write the plots and summary.tsv")
    args = parser.parse_args()

    strains = read_strains(args.strains)
    os.makedirs(args.out_dir, exist_ok=True)

    rows = []
    reason_rows = []
    severity_before_all, severity_after_all = {}, {}
    source_before_all, source_after_all = {}, {}

    for s in strains:
        original_path = os.path.join(args.original_dir, f"{s}_final_annotation.gff")
        filtered_path = os.path.join(args.filtered_dir, f"{s}_final_annotation_filtered.gff")
        report_path = os.path.join(args.filtered_dir, f"{s}_filter_report.tsv")

        for p in (original_path, filtered_path, report_path):
            if not os.path.isfile(p):
                sys.exit(f"ERROR: missing {p}")

        print(f"[{s}] loading...", file=sys.stderr)
        genes_before = parse_gff_for_qc(original_path, "protein_coding_gene")
        genes_after = parse_gff_for_qc(filtered_path, "gene")

        n_before = len(genes_before)
        n_after = len(genes_after)
        n_extra_before = sum(1 for g in genes_before.values() if g["extra_copy_number"])
        n_extra_after = sum(1 for g in genes_after.values() if g["extra_copy_number"])

        with open(report_path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                reason_rows.append((s, categorize_reason(row["reason"])))

        sev_before = severity_counts(genes_before, args.min_overlap_bp)
        sev_after = severity_counts(genes_after, args.min_overlap_bp)
        severity_before_all[s] = sev_before
        severity_after_all[s] = sev_after

        source_before_all[s] = source_counts(genes_before, strains[0])
        source_after_all[s] = source_counts(genes_after, strains[0])

        rows.append(
            {
                "strain": s,
                "n_genes_before": n_before,
                "n_genes_after": n_after,
                "n_removed": n_before - n_after,
                "n_extra_copy_before": n_extra_before,
                "n_extra_copy_after": n_extra_after,
                "n_overlap_pairs_before": sum(sev_before.values()),
                "n_overlap_pairs_after": sum(sev_after.values()),
                "n_contained_before": sev_before["contained"],
                "n_contained_after": sev_after["contained"],
            }
        )

    with open(os.path.join(args.out_dir, "summary.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # --- Plot 1: total gene counts ---
    fig, ax = plt.subplots(figsize=(max(8, len(strains) * 0.6), 5))
    x = range(len(strains))
    width = 0.35
    ax.bar([i - width / 2 for i in x], [r["n_genes_before"] for r in rows], width, label="before filtering")
    ax.bar([i + width / 2 for i in x], [r["n_genes_after"] for r in rows], width, label="after filtering")
    ax.set_xticks(list(x))
    ax.set_xticklabels(strains, rotation=45, ha="right")
    ax.set_ylabel("Genes")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "01_gene_counts.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 2: extra-copy genes ---
    fig, ax = plt.subplots(figsize=(max(8, len(strains) * 0.6), 5))
    ax.bar([i - width / 2 for i in x], [r["n_extra_copy_before"] for r in rows], width, label="before filtering")
    ax.bar([i + width / 2 for i in x], [r["n_extra_copy_after"] for r in rows], width, label="after filtering")
    ax.set_xticks(list(x))
    ax.set_xticklabels(strains, rotation=45, ha="right")
    ax.set_ylabel("Extra-copy genes")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "02_extra_copies.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 3: removal reasons, stacked ---
    categories = [label for _, label in REASON_CATEGORIES] + ["other"]
    reason_by_strain = {s: {c: 0 for c in categories} for s in strains}
    for s, cat in reason_rows:
        reason_by_strain[s][cat] += 1

    fig, ax = plt.subplots(figsize=(max(8, len(strains) * 0.6), 6))
    bottoms = [0] * len(strains)
    for cat in categories:
        values = [reason_by_strain[s][cat] for s in strains]
        if sum(values) == 0:
            continue
        ax.bar(strains, values, bottom=bottoms, label=cat)
        bottoms = [b + v for b, v in zip(bottoms, values)]
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("Genes removed")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.3), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "03_removal_reasons.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 4: overlap severity, before vs. after ---
    fig, axes = plt.subplots(2, 1, figsize=(max(8, len(strains) * 0.6), 9), sharex=True)
    sev_colors = {"minor": "#B0B0B0", "substantial": "#F2B035", "contained": "#C0392B"}
    for ax, data, title in [
        (axes[0], severity_before_all, "before filtering"),
        (axes[1], severity_after_all, "after filtering"),
    ]:
        bottoms = [0] * len(strains)
        for sev in ("minor", "substantial", "contained"):
            values = [data[s][sev] for s in strains]
            ax.bar(strains, values, bottom=bottoms, label=sev, color=sev_colors[sev])
            bottoms = [b + v for b, v in zip(bottoms, values)]
        ax.set_title(title)
        ax.set_ylabel("CDS-overlapping pairs")
    plt.setp(axes[1].get_xticklabels(), rotation=45, ha="right")
    axes[0].legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "04_overlap_severity.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 5: source composition, before vs. after ---
    src_labels = ["main", "original", "reference_rescued", "good_others", "rescued_others"]
    src_colors = {
        "main": "#F8766D",
        "original": "#00BA38",
        "reference_rescued": "#619CFF",
        "good_others": "#B379E0",
        "rescued_others": "#7A5195",
    }
    fig, axes = plt.subplots(2, 1, figsize=(max(8, len(strains) * 0.6), 9), sharex=True)
    for ax, data, title in [
        (axes[0], source_before_all, "before filtering"),
        (axes[1], source_after_all, "after filtering"),
    ]:
        bottoms = [0] * len(strains)
        for src in src_labels:
            values = [data[s][src] for s in strains]
            ax.bar(strains, values, bottom=bottoms, label=src, color=src_colors[src])
            bottoms = [b + v for b, v in zip(bottoms, values)]
        ax.set_title(title)
        ax.set_ylabel("Genes")
    plt.setp(axes[1].get_xticklabels(), rotation=45, ha="right")
    axes[0].legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.25), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "05_source_composition.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"\nDone - 5 plots + summary.tsv written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
