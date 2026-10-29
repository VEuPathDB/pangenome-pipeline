#!/usr/bin/env python3
"""
How similar are the models of one gene across genomes? For every basket with a reference gene R, and every genome T
with a native member in that basket, compare R's transfer onto T with T's own native model, in T's coordinates:
  identical          same introns (splice sites) and same CDS start and end
  same introns       same introns, CDS start or end differ
  different          introns differ (exon gained or lost, other splice sites), or T's gene is in pieces
and the in-frame coding overlap: bases coded by both in the same reading frame / bases coded by either (0-1).
Usage: python3 compare_models.py --baskets baskets.tsv --native-dir gff_originals --transfer-dir polished_nonredundant \
           --genomes genome_order.txt --out compare.tsv
"""
import argparse
import csv
import glob
import os
import re
from collections import Counter, defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def read_cds(path, want, primary_only=False):
    """gene ID -> (contig, strand, sorted CDS segments of its longest transcript)."""
    genes, tx, cds = {}, {}, defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or line.startswith("#"):
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if i and not p:
            if i.group(1) in want and not (primary_only and "extra_copy_number=" in f[8]
                                           and "extra_copy_number=0" not in f[8]):
                genes[i.group(1)] = (f[0], f[6])
        elif i and p and f[2] in ("mRNA", "transcript"):
            tx[i.group(1)] = p.group(1).split(",")[0]
        elif f[2] == "CDS" and p:
            cds[p.group(1).split(",")[0]].append((int(f[3]), int(f[4])))
    out = {}
    for t, segs in cds.items():
        g = tx.get(t, t)
        if g in genes and (g not in out or sum(b - a for a, b in segs) > sum(b - a for a, b in out[g][2])):
            out[g] = (genes[g][0], genes[g][1], sorted(segs))
    return out


def coded(segs, strand):
    """position -> frame (0,1,2) of every coding base, counted from the CDS start."""
    pos, n = {}, 0
    order = segs if strand == "+" else sorted(segs, reverse=True)
    for a, b in order:
        rng = range(a, b + 1) if strand == "+" else range(b, a - 1, -1)
        for x in rng:
            pos[x] = n % 3
            n += 1
    return pos


def introns(segs):
    return {(segs[i][1], segs[i + 1][0]) for i in range(len(segs) - 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baskets", required=True)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--transfer-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--reference", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    ref = a.reference or genomes[0]
    pref = {}
    bas = []
    for r in csv.DictReader(open(a.baskets), delimiter="\t"):
        mem = r["natives"].split(",")
        bas.append((r["basket"], mem))
    natives_all = {m for _, mem in bas for m in mem}
    native, owner = {}, {}
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))
        d = read_cds(path, natives_all)
        native.update(d)
        owner.update({k: g for k in d})
    refs = {m for _, mem in bas for m in mem if owner.get(m) == ref}
    trans = {}
    for g in genomes:
        if g == ref:
            continue
        p = os.path.join(a.transfer_dir, f"{ref}_onto_{g}_liftoff.gff3_polished")
        for k, v in read_cds(p, refs, primary_only=True).items():
            trans[(k, g)] = v
    rows, stats, scores = [], Counter(), defaultdict(list)
    for b, mem in bas:
        R = next((m for m in mem if owner.get(m) == ref), None)
        if not R:
            continue
        by = defaultdict(list)
        for m in mem:
            if owner.get(m) and owner[m] != ref:
                by[owner[m]].append(m)
        for t, ms in by.items():
            if (R, t) not in trans:
                continue
            c, st, rs = trans[(R, t)]
            if len(ms) > 1:
                cat = "different (in pieces)"
                own = sorted(s for m in ms for s in native[m][2])
            else:
                own = native[ms[0]][2]
                if native[ms[0]][0] != c or native[ms[0]][1] != st:
                    continue
                if introns(own) != introns(rs):
                    cat = "different"
                elif own[0][0] == rs[0][0] and own[-1][1] == rs[-1][1]:
                    cat = "identical"
                else:
                    cat = "same introns"
            p1, p2 = coded(rs, st), coded(own, st)
            same = sum(1 for x, f in p1.items() if p2.get(x) == f)
            score = same / len(set(p1) | set(p2)) if p1 else 0
            stats[cat] += 1
            scores[cat].append(score)
            rows.append([b, R, t, ",".join(ms), cat, f"{score:.3f}"])
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["basket", "reference_gene", "genome", "own_model", "category", "inframe_overlap"])
        w.writerows(rows)
    tot = sum(stats.values())
    print(f"reference model (its transfer) vs each genome's own model, {tot} comparisons:")
    for k in ("identical", "same introns", "different", "different (in pieces)"):
        s = sorted(scores[k])
        med = s[len(s) // 2] if s else 0
        print(f"  {k:24} {stats[k]:7} ({100 * stats[k] / tot:5.1f}%)   median in-frame overlap {med:.3f}")
    per = defaultdict(Counter)
    for b, R, t, ms, cat, sc in rows:
        per[t][cat] += 1
    print("\nper genome: % identical / same introns / different / in pieces")
    for t in genomes:
        if t in per:
            n = sum(per[t].values())
            print(f"  {t:10} " + "  ".join(f"{100 * per[t][k] / n:5.1f}" for k in
                                           ("identical", "same introns", "different", "different (in pieces)")))


if __name__ == "__main__":
    main()
