#!/usr/bin/env python3
"""Draw the native genes of chosen baskets (from build_baskets.py): one row per genome, CDS as boxes, drawn in the
genome's own coordinates from the leftmost member there, 5'->3' of the first member. Members in colour (one colour per
gene of the genome with 2+ members, else blue); partial members (pieces) in grey below.
Usage: python3 plot_baskets.py --baskets baskets.tsv --native-dir gff_originals --genomes genome_order.txt
                               (--ids B1,B2 | --random N --seed S) --out-prefix out"""
import argparse
import csv
import glob
import os
import random
import re
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
DESC_RE = re.compile(r"(?:^|;)description=([^;]+)")
COLS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2", "#937860"]


def genes_of(path, want):
    genes, tx, cds = {}, {}, defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or line.startswith("#"):
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if i and not p and i.group(1) in want:
            d = DESC_RE.search(f[8])
            genes[i.group(1)] = dict(c=f[0], st=f[6], S=int(f[3]), E=int(f[4]),
                                     desc=(d.group(1) if d else "").replace("%2C", ","))
        elif i and p and f[2] in ("mRNA", "transcript"):
            tx[i.group(1)] = p.group(1).split(",")[0]
        elif f[2] == "CDS" and p:
            cds[p.group(1).split(",")[0]].append((int(f[3]), int(f[4])))
    best = {}
    for t, segs in cds.items():
        g = tx.get(t, t)
        if g in genes and (g not in best or sum(b - a for a, b in segs) > sum(b - a for a, b in best[g])):
            best[g] = sorted(segs)
    for g in genes:
        genes[g]["cds"] = best.get(g, [])
    return genes


def draw(ax, y, segs, x0, strand, ref_end, color, h=0.32, alpha=1.0, gene_strand=None, tip=0.0):
    X = (lambda v: v - x0) if strand == "+" else (lambda v: ref_end - v)
    xs = sorted((min(X(a), X(b)), max(X(a), X(b))) for a, b in segs)
    if not xs:
        return
    ax.plot([xs[0][0], xs[-1][1]], [y, y], color=color, lw=1, alpha=alpha, zorder=1)
    for a, b in xs:
        ax.add_patch(Rectangle((a, y - h / 2), b - a, h, color=color, alpha=alpha, lw=0, zorder=2))
    if gene_strand and tip:   # arrowhead at the 3' end: right if on the drawn strand, left otherwise
        hh = max(h, 0.2)
        if gene_strand == strand:
            e = xs[-1][1]
            ax.add_patch(Polygon([(e, y - hh / 2), (e, y + hh / 2), (e + tip, y)], color=color, alpha=alpha, lw=0))
        else:
            e = xs[0][0]
            ax.add_patch(Polygon([(e, y - hh / 2), (e, y + hh / 2), (e - tip, y)], color=color, alpha=alpha, lw=0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baskets", required=True)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--ids", default="")
    ap.add_argument("--random", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    rows = {r["basket"]: r for r in csv.DictReader(open(a.baskets), delimiter="\t")}
    if a.ids:
        chosen = a.ids.split(",")
    else:
        pool = sorted(k for k, r in rows.items() if int(r["max_natives_one_genome"]) > 1)
        random.seed(a.seed)
        chosen = random.sample(pool, a.random)
    paths = {g: next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff"))))) for g in genomes}
    for b in chosen:
        r = rows[b]
        mem = r["natives"].split(",")
        pcs = [x for x in r["partial_members"].split(",") if x]
        info = {}
        for g in genomes:
            info.update({k: dict(v, g=g, piece=False) for k, v in genes_of(paths[g], set(mem)).items()})
            info.update({k: dict(v, g=g, piece=True) for k, v in genes_of(paths[g], set(pcs)).items()})
        rowsg = [g for g in genomes if any(v["g"] == g for v in info.values())]
        fig, ax = plt.subplots(figsize=(13, 0.62 * len(rowsg) + 1.6))
        xmax = 0
        span_all = max(max(v["E"] for v in info.values() if v["g"] == g) - min(v["S"] for v in info.values() if v["g"] == g)
                       for g in rowsg)
        tipw = span_all * 0.012
        for i, g in enumerate(rowsg):
            y = len(rowsg) - i
            here = sorted((k for k, v in info.items() if v["g"] == g), key=lambda k: info[k]["S"])
            members = [k for k in here if not info[k]["piece"]]
            first = info[members[0]] if members else info[here[0]]
            st = first["st"]
            x0 = min(info[k]["S"] for k in here)
            end = max(info[k]["E"] for k in here)
            multi = len(members) > 1
            ax.text(-150, y, g, ha="right", va="center", fontsize=10, fontweight="bold" if multi else "normal")
            lab = []
            for j, k in enumerate(members):
                col = COLS[j % len(COLS)] if multi else COLS[0]
                off = (j - (len(members) - 1) / 2) * 0.17 if multi else 0
                draw(ax, y + off, info[k]["cds"], x0, st, end, col, h=0.15 if multi else 0.3,
                     gene_strand=info[k]["st"], tip=tipw)
                lab.append((k, col, info[k]["st"]))
            for k in (k for k in here if info[k]["piece"]):
                draw(ax, y - 0.33, info[k]["cds"], x0, st, end, "#9AA5AE", h=0.1, alpha=0.8,
                     gene_strand=info[k]["st"], tip=tipw)
            xmax = max(xmax, end - x0)
            txt = "  ".join(f"{k}({s})" for k, _, s in lab)
            ax.text(1.02, y, txt, transform=ax.get_yaxis_transform(), fontsize=8, va="center", color="#52606D")
        desc = next((v["desc"] for k, v in info.items() if v["g"] == genomes[0] and not v["piece"]), "") or \
            next((v["desc"] for v in info.values() if v["desc"]), "")
        ax.set_title(f"Basket {b}: {desc}   (bold: genome with 2+ member genes; grey below: pieces; arrow = 3' end)", loc="left",
                     fontsize=11)
        ax.set_xlim(-200, xmax * 1.03)
        ax.set_ylim(0.3, len(rowsg) + 0.7)
        ax.set_yticks([])
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("bp from the leftmost member gene, 5'->3' of the first member")
        out = f"{a.out_prefix}_{b}.png"
        fig.savefig(out, dpi=160, bbox_inches="tight")
        plt.close(fig)
        print(out)


if __name__ == "__main__":
    main()
