#!/usr/bin/env python3
"""Draw one panel gene and the panel genes it overlaps, in every genome: each gene's model there (its own gene in its
genome, its primary Liftoff transfer elsewhere, coverage >= 0.8), CDS as boxes, arrow at the 3' end. Rows where the
genes overlap are marked.
Usage: python3 plot_overlap_case.py --baskets baskets_v2.tsv --native-dir gff_originals --transfer-dir polished_nonredundant
           --domain-dir annotation_check_with_interpro --genomes genome_order.txt --gene ID --out fig.png"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
COV_RE = re.compile(r"(?:^|;)coverage=([\d.]+)")
DESC_RE = re.compile(r"(?:^|;)description=([^;]+)")
COLS = ["#4C72B0", "#DD8452", "#55A868", "#8172B2"]


def models(path, want, primary=False):
    genes, tx, cds = {}, {}, defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or line.startswith("#"):
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if i and not p and i.group(1) in want:
            if primary and "extra_copy_number=" in f[8] and "extra_copy_number=0" not in f[8]:
                continue
            c = COV_RE.search(f[8])
            d = DESC_RE.search(f[8])
            genes[i.group(1)] = (f[0], f[6], float(c.group(1)) if c else 1.0, (d.group(1) if d else "").replace("%2C", ","))
        elif i and p and f[2] in ("mRNA", "transcript"):
            tx[i.group(1)] = p.group(1).split(",")[0]
        elif f[2] == "CDS" and p:
            cds[p.group(1).split(",")[0]].append((int(f[3]), int(f[4])))
    out = {}
    for t, segs in cds.items():
        g = tx.get(t, t)
        if g in genes and (g not in out or sum(b - a for a, b in segs) > sum(b - a for a, b in out[g][4])):
            out[g] = genes[g] + (sorted(segs),)
    return out


def main():
    ap = argparse.ArgumentParser()
    for k in ("baskets", "native-dir", "transfer-dir", "domain-dir", "genomes", "gene", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip()]
    winners = {w for r in csv.DictReader(open(a.baskets), delimiter="\t") for w in r["winners"].split(";")}
    owner = {}
    nat = {}
    for g in genomes:
        d = models(next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff"))))), winners)
        for k, v in d.items():
            owner[k], nat[k] = g, v
    dom = {}
    for g in genomes:
        for r in csv.DictReader(open(os.path.join(a.domain_dir, f"{g}_gene_classification.tsv")), delimiter="\t"):
            dom[r["gene_id"]] = r["has_domain"] == "True"
    per = defaultdict(dict)   # genome -> winner -> model
    for w, v in nat.items():
        per[owner[w]][w] = v
    for path in sorted(glob.glob(os.path.join(a.transfer_dir, "*_onto_*"))):
        src, tgt = os.path.basename(path).split("_liftoff")[0].split("_onto_")
        if tgt not in genomes:
            continue
        want = {w for w in winners if owner.get(w) == src}
        for k, v in models(path, want, primary=True).items():
            if v[2] >= 0.8:
                per[tgt][k] = v
    G = a.gene
    partners = set()
    for t in genomes:
        m = per[t].get(G)
        if not m:
            continue
        s, e = m[4][0][0], m[4][-1][1]
        for k, v in per[t].items():
            if k != G and v[0] == m[0] and min(e, v[4][-1][1]) - max(s, v[4][0][0]) + 1 >= 20:
                partners.add(k)
    genes = [G] + sorted(partners)
    fig, ax = plt.subplots(figsize=(13, 0.62 * len(genomes) + 2))
    for i, t in enumerate(genomes):
        y = len(genomes) - i
        ms = [(k, per[t].get(k)) for k in genes]
        ms = [(k, m) for k, m in ms if m]
        ax.text(-60, y, t, ha="right", va="center", fontsize=10)
        if not ms:
            continue
        x0 = min(m[4][0][0] for _, m in ms)
        x1 = max(m[4][-1][1] for _, m in ms)
        st = (per[t].get(G) or ms[0][1])[1]
        X = (lambda v: v - x0) if st == "+" else (lambda v: x1 - v)
        mg = per[t].get(G)
        over = False
        for j, k in enumerate(genes):
            m = per[t].get(k)
            if not m:
                continue
            yy = y + (0.17 if j == 0 else -0.17)
            col = COLS[j % len(COLS)]
            xs = sorted((min(X(p), X(q)), max(X(p), X(q))) for p, q in m[4])
            ax.plot([xs[0][0], xs[-1][1]], [yy, yy], color=col, lw=1)
            for p, q in xs:
                ax.add_patch(Rectangle((p, yy - 0.12), q - p, 0.24, color=col, lw=0))
            tip = (x1 - x0) * 0.015 + 15
            if m[1] == st:
                ax.add_patch(Polygon([(xs[-1][1], yy - 0.15), (xs[-1][1], yy + 0.15), (xs[-1][1] + tip, yy)], color=col, lw=0))
            else:
                ax.add_patch(Polygon([(xs[0][0], yy - 0.15), (xs[0][0], yy + 0.15), (xs[0][0] - tip, yy)], color=col, lw=0))
            if j and mg and m[0] == mg[0] and min(mg[4][-1][1], m[4][-1][1]) - max(mg[4][0][0], m[4][0][0]) + 1 >= 20:
                over = True
        note = "overlap" if over else ("" if mg else f"{G} absent")
        ax.text(1.01, y, note, transform=ax.get_yaxis_transform(), fontsize=9, va="center",
                color="#D62728" if over else "#52606D", fontweight="bold" if over else "normal")
    title = []
    for j, k in enumerate(genes):
        desc = nat[k][3][:40]
        title.append(f"{'blue' if j == 0 else ['orange', 'green', 'purple'][(j - 1) % 3]}: {k} ({owner[k]}; "
                     f"{'domain' if dom.get(k) else 'no domain'}; {desc})")
    ax.set_title("\n".join(title) + "\narrow = 3' end", loc="left", fontsize=9.5)
    ax.set_yticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_ylim(0.3, len(genomes) + 0.7)
    ax.set_xlabel("bp in each genome, from the leftmost of the genes drawn")
    fig.savefig(a.out, dpi=160, bbox_inches="tight")
    print(a.out, genes)


if __name__ == "__main__":
    main()
