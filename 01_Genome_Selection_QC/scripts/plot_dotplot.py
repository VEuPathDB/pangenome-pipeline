#!/usr/bin/env python3
"""
Synteny dot-plot from a minimap2 PAF: reference on X, query on Y, with the query
contigs ordered and oriented by where they land on the reference so collinear
pieces stack along the diagonal.

    plot_dotplot.py aln.paf --ref-name Af293 --query-name GCA_024220425.1 -o out.png

Blue segment  = forward-strand alignment
Red segment   = reverse strand (inversion)
Only alignments >= --min-aln bp are drawn, so large-scale rearrangements stand out
instead of drowning in short hits.
"""
import argparse
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_paf(path, min_aln):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            aln = int(f[10])
            if aln < min_aln:
                continue
            rows.append(dict(
                q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
                t=f[5], tlen=int(f[6]), ts=int(f[7]), te=int(f[8])))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf")
    ap.add_argument("--ref-name", default="reference")
    ap.add_argument("--query-name", default="query")
    ap.add_argument("--min-aln", type=int, default=2000)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    rows = read_paf(a.paf, a.min_aln)
    if not rows:
        sys.exit(f"no alignments >= {a.min_aln} bp in {a.paf}")

    # X axis: reference sequences, biggest first, laid end to end
    tlen = {r["t"]: r["tlen"] for r in rows}
    t_order = sorted(tlen, key=lambda t: -tlen[t])
    t_off, cum = {}, 0
    for t in t_order:
        t_off[t] = cum
        cum += tlen[t]
    x_total = cum
    t_rank = {t: i for i, t in enumerate(t_order)}

    # Y axis: query contigs, placed after the reference sequence they align to
    # most, then by mean position on it
    qlen = {r["q"]: r["qlen"] for r in rows}
    aligned = defaultdict(lambda: defaultdict(float))   # q -> t -> bp
    pos_sum = defaultdict(lambda: defaultdict(float))
    pos_n = defaultdict(lambda: defaultdict(int))
    for r in rows:
        span = r["te"] - r["ts"]
        aligned[r["q"]][r["t"]] += span
        pos_sum[r["q"]][r["t"]] += (r["ts"] + r["te"]) / 2
        pos_n[r["q"]][r["t"]] += 1

    def key(q):
        main_t = max(aligned[q], key=aligned[q].get)
        mean_pos = pos_sum[q][main_t] / max(1, pos_n[q][main_t])
        return (t_rank[main_t], mean_pos)

    q_order = sorted(qlen, key=key)
    q_off, cum = {}, 0
    for q in q_order:
        q_off[q] = cum
        cum += qlen[q]
    y_total = cum

    ratio = min(2.0, max(0.35, y_total / x_total))
    fig, ax = plt.subplots(figsize=(15, 15 * ratio))
    for r in rows:
        x0 = t_off[r["t"]] + r["ts"]
        x1 = t_off[r["t"]] + r["te"]
        if r["strand"] == "+":
            y0, y1, c = q_off[r["q"]] + r["qs"], q_off[r["q"]] + r["qe"], "#3b5bdb"
        else:
            y0, y1, c = q_off[r["q"]] + r["qe"], q_off[r["q"]] + r["qs"], "#e03131"
        ax.plot([x0, x1], [y0, y1], "-", color=c, lw=1.1, alpha=0.8)

    for t in t_order:
        ax.axvline(t_off[t], color="0.75", lw=0.6)
    for q in q_order:
        ax.axhline(q_off[q], color="0.9", lw=0.4)

    ax.set_xlim(0, x_total)
    ax.set_ylim(0, y_total)
    ax.set_xticks([t_off[t] + tlen[t] / 2 for t in t_order])
    ax.set_xticklabels(t_order, rotation=90, fontsize=7)
    big_q = [q for q in q_order if qlen[q] > y_total * 0.01]
    ax.set_yticks([q_off[q] + qlen[q] / 2 for q in big_q])
    ax.set_yticklabels(big_q, fontsize=7)
    ax.set_xlabel(a.ref_name)
    ax.set_ylabel(a.query_name)
    ax.set_title(f"{a.query_name}  vs  {a.ref_name}     "
                 f"(blue = same strand, red = inversion; alns >= {a.min_aln} bp)")
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f"wrote {a.out}  ({len(rows)} alignments, {len(q_order)} query contigs)")


if __name__ == "__main__":
    main()
