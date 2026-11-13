#!/usr/bin/env python3
"""
Plots of the basket annotation (basket_annotation_job.sh output), one bar per genome, the modes side by side.

  <label>_model_origin.png   where each gene's model comes from: the genome's own model (kept, or kept because no
                             transfer fitted), the panel model (replaced, filled, alternative), repaired by miniprot,
                             or left with an invalid ORF
  <label>_orf_quality.png    the final genes: valid ORF with / without a protein domain, and the genes still without a
                             valid ORF, by what breaks them (miniprot's diagnosis)
  <label>_gene_counts.png    genes in each genome's own annotation and in the final GFF of each mode
  <label>_panel_genes.png    the panel genes: kept (with / without a reference gene) and removed (overlap, curation)
  <label>_miniprot.png       the miniprot results, all genomes together, per mode
  <label>_summary.tsv        the numbers behind the plots

Usage: python3 plot_basket_annotation.py --result-dir basket_annotation --genomes genome_order.txt \\
           --native-dir gff_for_liftoff --label Toxo --out-dir basket_annotation/plots [--modes "keep_similar prefer_panel_model"]
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_baskets import read_gff   # noqa: E402

ORIGINS = [("original", "own model kept", "#55A868"),
           ("original_divergent", "own model kept (no transfer fitted)", "#A6D8A8"),
           ("original_other_place", "own gene kept at another place of its basket", "#2E7D32"),
           ("replaced", "replaced by the panel model", "#4C72B0"),
           ("filled", "filled with the panel model (genome had no gene)", "#8172B2"),
           ("alternative", "another model, closest to the panel model", "#CCB974"),
           ("rescued", "rescued: a partial valid model (>= 50% of the panel model)", "#DD8452"),
           ("repaired_miniprot", "repaired by miniprot", "#64B5CD"),
           ("target_invalid_orf", "panel model, invalid ORF", "#C44E52")]
MODE_TITLE = {"keep_similar": "keep_similar: own model kept when similar to the panel model",
              "prefer_panel_model": "prefer_panel_model: the valid model closest to the panel model"}
QUALITY = [("valid_domain", "valid ORF, domain", "#2C7FB8"),
           ("valid_nodomain", "valid ORF, no domain", "#A6CEE3"),
           ("repaired", "repaired by miniprot", "#41AB5D"),
           ("assembly", "N in the CDS or near a contig end", "#BDBDBD"),
           ("frameshift", "frameshift", "#E31A1C"),
           ("stop_codon", "internal stop", "#FB6A4A"),
           ("start_or_stop_lost", "start or stop codon lost", "#FC9272"),
           ("noncanonical", "repair needs a non-canonical intron", "#FDD0A2"),
           ("partial", "partial / no miniprot hit", "#FDBB84"),
           ("selenoprotein", "selenoprotein", "#BCBDDC"),
           ("pseudogene", "pseudogene (annotated)", "#6A51A3"),
           ("other_invalid", "invalid ORF, other", "#969696")]


def gff_genes(path):
    out = []
    for line in open(path):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "gene":
            out.append(f[8])
    return out


def quality(attrs):
    if "model_origin=repaired_miniprot" in attrs:
        return "repaired"
    if "orf_invalid=True" not in attrs:
        return "valid_domain" if "has_domain=True" in attrs else "valid_nodomain"
    if "selenoprotein=True" in attrs:
        return "selenoprotein"
    if "pseudogene=True" in attrs:
        return "pseudogene"
    if "frameshift=" in attrs:
        return "frameshift"
    if "stop_codon=" in attrs:
        return "stop_codon"
    if "start_codon_lost=True" in attrs or "stop_codon_lost=True" in attrs:
        return "start_or_stop_lost"
    if "orf_status=invalid_assembly_gap" in attrs or "orf_status=invalid_near_contig_end" in attrs:
        return "assembly"
    if "miniprot_repair_noncanonical_intron=True" in attrs:
        return "noncanonical"
    if "partial_alignment=True" in attrs or "no_miniprot_hit=True" in attrs:
        return "partial"
    return "other_invalid"


def stacked(ax, labels, rows, keys, title):
    left = [0] * len(labels)
    for k, lab, col in keys:
        v = [r.get(k, 0) for r in rows]
        ax.barh(labels, v, left=left, color=col, label=lab, edgecolor="white", linewidth=0.4)
        left = [a + b for a, b in zip(left, v)]
    ax.set_title(title, loc="left", fontsize=11)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def save(fig, base):
    fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", base + ".png")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--result-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--modes", default="prefer_panel_model")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    modes = a.modes.split()
    base = os.path.join(a.out_dir, a.label)

    def gdir(m):
        """the final GFFs of mode m: after the overlap rule, else after miniprot, else the build."""
        for d in (f"gff_{m}_final", f"gff_{m}_miniprot", f"gff_{m}"):
            if os.path.isdir(os.path.join(a.result_dir, d)):
                return os.path.join(a.result_dir, d)
        return os.path.join(a.result_dir, f"gff_{m}")

    origin, qual, count = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    nonref = defaultdict(dict)
    for m in modes:
        for g in genomes:
            attrs = gff_genes(os.path.join(gdir(m), f"{g}_basket_annotation.gff"))
            origin[m][g] = Counter(re.search(r"model_origin=([^;]+)", x).group(1) for x in attrs)
            qual[m][g] = Counter(quality(x) for x in attrs)
            count[m][g] = len(attrs)
            nonref[m][g] = sum("not_in_reference_annotation=True" in x for x in attrs)
    own = {g: len(read_gff(next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff"))))))) for g in genomes}

    h = 0.42 * len(genomes) + 2.2
    # 1. model origin
    fig, axes = plt.subplots(1, len(modes), figsize=(7 * len(modes), h), sharey=True, squeeze=False)
    for ax, m in zip(axes[0], modes):
        stacked(ax, genomes, [origin[m][g] for g in genomes], ORIGINS, MODE_TITLE.get(m, m))
        ax.set_xlabel("genes")
    axes[0][0].invert_yaxis()
    h_, l_ = axes[0][0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    save(fig, base + "_model_origin")
    # 2. ORF quality (in %)
    fig, axes = plt.subplots(1, len(modes), figsize=(7 * len(modes), h), sharey=True, squeeze=False)
    for ax, m in zip(axes[0], modes):
        rows = [{k: 100 * v / max(count[m][g], 1) for k, v in qual[m][g].items()} for g in genomes]
        stacked(ax, genomes, rows, QUALITY, MODE_TITLE.get(m, m))
        ax.set_xlim(0, 100)
        ax.set_xlabel("% of the genes")
    axes[0][0].invert_yaxis()
    h_, l_ = axes[0][0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    save(fig, base + "_orf_quality")
    # 3. gene counts
    fig, ax = plt.subplots(figsize=(9, h))
    n = len(modes) + 1
    w = 0.8 / n
    ys = range(len(genomes))
    ax.barh([y - 0.4 + w / 2 for y in ys], [own[g] for g in genomes], height=w, color="#BBBBBB", label="own annotation")
    cols = ["#4C72B0", "#DD8452", "#55A868"]
    for i, m in enumerate(modes):
        ax.barh([y - 0.4 + w / 2 + w * (i + 1) for y in ys], [count[m][g] for g in genomes], height=w,
                color=cols[i % 3], label=f"final, {m}")
    ax.set_yticks(list(ys))
    ax.set_yticklabels(genomes)
    ax.invert_yaxis()
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel("protein-coding genes")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=n, frameon=False)
    save(fig, base + "_gene_counts")
    # 4. panel genes
    winners = list(csv.DictReader(open(os.path.join(a.result_dir, "winners.tsv")), delimiter="\t"))
    removed = {r["panel_gene"]: r["reason"] for r in
               csv.DictReader(open(os.path.join(gdir(modes[0]), "removed_panel_genes.tsv")), delimiter="\t")}
    pc = Counter()
    for r in winners:
        ref = r["is_reference"] == "1"
        why = removed.get(r["gene"])
        if why:
            pc["removed: " + ("overlap" if why.startswith("overlap") else "removed by the reference curation")] += 1
        else:
            pc["kept: reference gene" if ref else "kept: not in the reference annotation"] += 1
    keys = ["kept: reference gene", "kept: not in the reference annotation", "removed: overlap",
            "removed: removed by the reference curation"]
    fig, ax = plt.subplots(figsize=(9, 2.8))
    vals = [pc[k] for k in keys]
    colors = ["#4C72B0", "#8172B2", "#C44E52", "#937860"]
    ax.barh(keys, vals, color=colors)
    for y, v in enumerate(vals):
        ax.text(v + max(vals) * 0.01, y, f"{v:,}", va="center", fontsize=9)
    ax.invert_yaxis()
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_xlabel("panel genes")
    save(fig, base + "_panel_genes")
    # 5. miniprot
    res = {}
    for m in modes:
        p = os.path.join(a.result_dir, f"gff_{m}_miniprot", "miniprot_results.tsv")
        if os.path.exists(p):
            res[m] = Counter(r["result"] for r in csv.DictReader(open(p), delimiter="\t"))
    if res:
        cats = sorted({k for c in res.values() for k in c}, key=lambda k: -sum(c.get(k, 0) for c in res.values()))
        fig, ax = plt.subplots(figsize=(9, 0.45 * len(cats) + 1.5))
        w = 0.8 / len(res)
        for i, (m, c) in enumerate(res.items()):
            ax.barh([y - 0.4 + w / 2 + w * i for y in range(len(cats))], [c.get(k, 0) for k in cats], height=w,
                    color=cols[i % 3], label=m)
        ax.set_yticks(range(len(cats)))
        ax.set_yticklabels([k.replace("_", " ") for k in cats])
        ax.invert_yaxis()
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("genes with an invalid ORF, all genomes")
        ax.legend(frameon=False)
        save(fig, base + "_miniprot")
    with open(base + "_summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        okeys = [k for k, _, _ in ORIGINS]
        qkeys = [k for k, _, _ in QUALITY]
        w.writerow(["mode", "genome", "own_annotation", "final_genes", "not_in_reference_annotation"] + okeys + qkeys)
        for m in modes:
            for g in genomes:
                w.writerow([m, g, own[g], count[m][g], nonref[m][g]] + [origin[m][g].get(k, 0) for k in okeys]
                           + [qual[m][g].get(k, 0) for k in qkeys])
    print("wrote", base + "_summary.tsv")


if __name__ == "__main__":
    main()
