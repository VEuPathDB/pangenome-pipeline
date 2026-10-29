#!/usr/bin/env python3
"""
Why genes did not become panel genes in step 05: the genes in not_eligible.tsv (written by baskets_v2.py), by reason
and by the genome they come from. These genes are not in the final GFFs unless another member of their basket won.

Reasons (the filter of the winner choice):
  PSAURON < 0.9, sequence in the reference   a gene the reference does not annotate although its sequence is there
                                             (the stricter cut-off)
  PSAURON < 0.7                              any other gene
  protein shorter than 100 aa
  no valid ORF in any genome

Outputs (in --out-dir): <label>_not_eligible.png and <label>_not_eligible.tsv (genome x reason counts)

Usage:
    python3 plot_not_eligible.py --result-dir basket_annotation --label Toxo --out-dir basket_annotation/plots
"""
import argparse
import csv
import os
import re
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REASONS = [("ref_seq_psauron", "PSAURON < 0.9 (not annotated in the reference, sequence there)", "#2a6f97"),
           ("psauron", "PSAURON < 0.7", "#8fb8d0"),
           ("short", "protein shorter than 100 aa", "#e9a03b"),
           ("no_valid_orf", "no valid ORF in any genome", "#c0392b"),
           ("other", "other", "#9a9a9a")]


def reason_key(r):
    if r.startswith("PSAURON") and "reference" in r:
        return "ref_seq_psauron"
    if r.startswith("PSAURON"):
        return "psauron"
    if r.startswith("protein shorter"):
        return "short"
    if r.startswith("no valid ORF"):
        return "no_valid_orf"
    return "other"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--result-dir", required=True, help="step 05 output (not_eligible.tsv)")
    ap.add_argument("--label", default="Toxo")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(os.path.join(a.result_dir, "not_eligible.tsv")), delimiter="\t"))
    c = Counter((r["genome"], reason_key(r["reason"])) for r in rows)
    tot = Counter(reason_key(r["reason"]) for r in rows)
    reasons = [x for x in REASONS if tot[x[0]]]
    genomes = sorted({r["genome"] for r in rows}, key=lambda g: sum(c[(g, k)] for k, _, _ in reasons))

    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, f"{a.label}_not_eligible.tsv"), "w") as fh:
        fh.write("genome\t" + "\t".join(k for k, _, _ in reasons) + "\ttotal\n")
        for g in genomes[::-1]:
            fh.write(g + "\t" + "\t".join(str(c[(g, k)]) for k, _, _ in reasons)
                     + f"\t{sum(c[(g, k)] for k, _, _ in reasons)}\n")
        fh.write("all\t" + "\t".join(str(tot[k]) for k, _, _ in reasons) + f"\t{len(rows)}\n")

    fig, (x0, x1) = plt.subplots(1, 2, figsize=(14, 0.38 * len(genomes) + 2.6),
                                 gridspec_kw=dict(width_ratios=[1, 2.2]))
    # left: all genomes together, one bar per reason
    ys = range(len(reasons))[::-1]
    x0.barh(list(ys), [tot[k] for k, _, _ in reasons], color=[col for _, _, col in reasons], height=0.6)
    for y, (k, _, _) in zip(ys, reasons):
        x0.text(tot[k], y, f" {tot[k]:,} ({100 * tot[k] / len(rows):.0f}%)", va="center", fontsize=9)
    x0.set_yticks(list(ys), [re.sub(r" \(", "\n(", lab) for _, lab, _ in reasons], fontsize=9)
    x0.set_xlabel("genes")
    x0.set_xlim(0, max(tot.values()) * 1.35)
    x0.set_title(f"All genomes: {len(rows):,} genes", loc="left", fontsize=11)
    # right: per genome, stacked
    left = [0] * len(genomes)
    for k, lab, col in reasons:
        v = [c[(g, k)] for g in genomes]
        x1.barh(genomes, v, left=left, color=col, height=0.7, label=lab, edgecolor="white", linewidth=0.8)
        left = [l + n for l, n in zip(left, v)]
    for i, n in enumerate(left):
        x1.text(n, i, f" {n:,}", va="center", fontsize=8)
    x1.set_xlabel("genes")
    x1.set_title("By genome (the genome the gene comes from)", loc="left", fontsize=11)
    x1.legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2)
    for x in (x0, x1):
        x.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"{a.label}: genes that could not become panel genes in step 05, by reason", x=0.01, ha="left",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(a.out_dir, f"{a.label}_not_eligible.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
