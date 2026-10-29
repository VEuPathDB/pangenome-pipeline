#!/usr/bin/env python3
"""
Overlaps between panel genes (winners of baskets_v2.py), and how many would be removed by step 07's pair rules at
different vote thresholds. A winner's model in genome t is its own gene (its genome) or its primary transfer with
coverage >= --min-cov. Two winners overlap in t when their CDS spans (first to last CDS, introns included) overlap
>= --min-bp, on either strand.
Pair rules (who "would drop" in one genome): both reference -> neither; reference vs no domain -> the other;
reference vs domain -> neither; only one with a domain -> the one without; both with a domain -> neither; neither with
a domain -> fewer genomes with a valid ORF, then fewer genomes, then the shorter CDS.
A gene is removed when it would drop in a share of the genomes where it exists >= the threshold.
Domains: annotation_check_with_interpro/<genome>_gene_classification.tsv (has_domain of each native gene).
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_baskets import read_gff   # noqa: E402

VO = re.compile(r"valid_ORFs?=(\d+)")
COV = re.compile(r"(?:^|;)coverage=([\d.]+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baskets", required=True)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--transfer-dir", required=True)
    ap.add_argument("--domain-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--reference", default="")
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--min-bp", type=int, default=20)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    ref = a.reference or genomes[0]
    winners = [w for r in csv.DictReader(open(a.baskets), delimiter="\t") for w in r["winners"].split(";")]
    wset = set(winners)
    owner, native = {}, {}
    for g in genomes:
        for k, v in read_gff(next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))).items():
            if k in wset:
                owner[k], native[k] = g, v
    dom = {}
    for g in genomes:
        for r in csv.DictReader(open(os.path.join(a.domain_dir, f"{g}_gene_classification.tsv")), delimiter="\t"):
            dom[r["gene_id"]] = r["has_domain"] == "True"
    models = defaultdict(list)   # genome -> (contig, strand, s, e, winner)
    nvalid, npres = Counter(), Counter()
    for w in winners:
        c, st, s, e, attr = native[w]
        models[owner[w]].append((c, st, s, e, w))
        npres[w] += 1
        nvalid[w] += 1
    for path in sorted(glob.glob(os.path.join(a.transfer_dir, "*_onto_*"))):
        src, tgt = os.path.basename(path).split("_liftoff")[0].split("_onto_")
        if tgt not in genomes:
            continue
        for k, (c, st, s, e, attr) in read_gff(path, primary_only=True).items():
            if k in wset and owner.get(k) == src:
                cv = COV.search(attr)
                if cv and float(cv.group(1)) >= a.min_cov:
                    models[tgt].append((c, st, s, e, k))
                    npres[k] += 1
                    m = VO.search(attr)
                    nvalid[k] += bool(m and m.group(1) != "0")
    pairs = defaultdict(set)
    for t, L in models.items():
        L.sort(key=lambda x: (x[0], x[2]))
        act = []
        for c, st, s, e, w in L:
            act = [x for x in act if x[0] == c and x[3] >= s]
            for c2, st2, s2, e2, w2 in act:
                if w2 != w and min(e, e2) - max(s, s2) + 1 >= a.min_bp:
                    pairs[tuple(sorted((w, w2)))].add(t)
            act.append((c, st, s, e, w))
    present_in = defaultdict(set)
    for t, L in models.items():
        for *_, w in L:
            present_in[w].add(t)

    def loser(x, y):
        rx, ry = owner[x] == ref, owner[y] == ref
        dx, dy = dom.get(x, False), dom.get(y, False)
        if rx and ry:
            return None
        if rx or ry:
            r_, o_ = (x, y) if rx else (y, x)
            return None if dom.get(o_, False) else o_
        if dx != dy:
            return y if dx else x
        if dx and dy:
            return None
        kx = (nvalid[x] / max(npres[x], 1), npres[x], native[x][3] - native[x][2])
        ky = (nvalid[y] / max(npres[y], 1), npres[y], native[y][3] - native[y][2])
        return x if kx < ky else y

    votes = defaultdict(set)   # gene -> genomes where it would drop
    pair_share = defaultdict(float)   # gene -> highest share over its pairs, counting only genomes with both genes
    stats = Counter()
    for (x, y), gs in pairs.items():
        lo = loser(x, y)
        stats["overlapping pairs"] += 1
        stats["pairs where nobody drops (reference pair, or domain protects)" if lo is None else "pairs with a loser"] += 1
        if lo:
            votes[lo] |= gs
            both = present_in[x] & present_in[y]
            pair_share[lo] = max(pair_share[lo], len(gs) / max(len(both), 1))
    rows = []
    for g, gs in votes.items():
        rows.append((g, owner[g], len(gs), len(present_in[g]), len(gs) / max(len(present_in[g]), 1),
                     round(pair_share.get(g, 0), 3)))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "genome", "would_drop_in", "present_in", "share", "share_where_both_exist"])
        w.writerows(sorted(rows, key=lambda r: -r[4]))
    print(f"panel genes: {len(winners)}")
    for k, v in stats.items():
        print(f"  {v:6}  {k}")
    print(f"genes that would drop in at least one genome: {len(rows)}")
    for lab, thr, strict in (("in any genome", 0.0, True), ("in >= 1/3 of its genomes", 1 / 3, False),
                             ("in more than half (step 07 rule)", 0.5, True), ("in all its genomes", 1.0, False)):
        sel = [r for r in rows if (r[4] > thr if strict else r[4] >= thr)]
        nonref = sum(1 for r in sel if r[1] != ref)
        print(f"  removed if it would drop {lab:34}: {len(sel):5}  (from other annotations: {nonref}, reference genes: {len(sel) - nonref})")
    print("\ncounting only the genomes where BOTH genes of the pair exist:")
    for lab, thr, strict in (("in any of them", 0.0, True), ("in >= 1/3 of them", 1 / 3, False),
                             ("in more than half", 0.5, True), ("in all of them", 1.0, False)):
        sel = [g for g, v in pair_share.items() if (v > thr if strict else v >= thr)]
        print(f"  removed if it would drop {lab:20}: {len(sel):5}  (reference genes: {sum(1 for g in sel if owner[g] == ref)})")
    ps = Counter(min(int(v * 10), 9) for v in pair_share.values())
    print("  share (0.0-0.1 ... 0.9-1.0):", [ps[i] for i in range(10)])
    shares = Counter(min(int(r[4] * 10), 9) for r in rows)
    print("share of genomes where the gene would drop (0.0-0.1 ... 0.9-1.0):",
          [shares[i] for i in range(10)])


if __name__ == "__main__":
    main()
