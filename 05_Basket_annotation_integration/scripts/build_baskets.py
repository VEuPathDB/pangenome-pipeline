#!/usr/bin/env python3
"""
Test of the gene "baskets" (loci of the panel) for the revised selection - steps 1-3 of REVISION_PLAN.md.

Nodes are the native protein-coding genes of every genome. Two native genes a (genome A) and b (genome B) are linked
when a's Liftoff transfer onto B covers b and b's transfer onto A covers a: CDS spans on the same strand, overlapping by
>= --min-frac of each. Linked genes form a basket (union-find). A native gene that lies >= --partial-frac inside a
member's transfer but is much shorter (< --min-frac of it) is listed as a partial member (a piece), and links nothing.
Only primary Liftoff copies (extra_copy_number=0) are used.

Usage:
    python3 build_baskets.py --native-dir gff_originals --transfer-dir polished_nonredundant \
        --genomes genome_order.txt --out-dir baskets_test
"""
import argparse
import bisect
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def read_gff(path, primary_only=False):
    """gene ID -> (contig, strand, cds_start, cds_end, attributes of the gene line), protein-coding genes with a CDS."""
    genes, tx2gene, cds = {}, {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t = f[2]
            if t in ("protein_coding_gene", "gene"):
                m = ID_RE.search(f[8])
                if m and "Parent=" not in f[8]:
                    if primary_only and "extra_copy_number=0" not in f[8] and "extra_copy_number=" in f[8]:
                        continue
                    genes[m.group(1)] = (f[0], f[6], f[8])
            elif t in ("mRNA", "transcript"):
                i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
                if i and p:
                    tx2gene[i.group(1)] = p.group(1).split(",")[0]
            elif t == "CDS":
                p = PARENT_RE.search(f[8])
                if p:
                    s, e = int(f[3]), int(f[4])
                    k = p.group(1).split(",")[0]
                    if k in cds:
                        a, b = cds[k]
                        cds[k] = (min(a, s), max(b, e))
                    else:
                        cds[k] = (s, e)
    out = {}
    for k, (s, e) in cds.items():
        g = tx2gene.get(k, k)
        if g in genes:
            c, st, attr = genes[g]
            if g in out:
                s, e = min(s, out[g][2]), max(e, out[g][3])
            out[g] = (c, st, s, e, attr)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--transfer-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--reference", default="")
    ap.add_argument("--min-frac", type=float, default=0.5)
    ap.add_argument("--partial-frac", type=float, default=0.8)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(args.genomes) if l.strip() and not l.startswith("#")]
    ref = args.reference or genomes[0]

    # natives
    native, idx = {}, {}
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(args.native_dir, g + ".gff3")) + glob.glob(os.path.join(args.native_dir, g + ".gff")))), None)
        genes = read_gff(path)
        native[g] = genes
        by = defaultdict(list)
        for gid, (c, st, s, e, _) in genes.items():
            by[(c, st)].append((s, e, gid))
        for k in by:
            by[k].sort()
        idx[g] = {k: ([x[0] for x in v], v) for k, v in by.items()}
        print(f"{g}: {len(genes)} native protein-coding genes", file=sys.stderr)
    owner = {gid: g for g in genomes for gid in native[g]}

    # transfers: directed hits a -> b, partial pieces
    hit = set()
    partial = defaultdict(set)   # basket member a -> pieces b in another genome
    n_tr = 0
    for path in sorted(glob.glob(os.path.join(args.transfer_dir, "*_onto_*"))):
        src, tgt = os.path.basename(path).split("_liftoff")[0].split("_onto_")
        if src not in native or tgt not in native:
            continue
        for gid, (c, st, s, e, _) in read_gff(path, primary_only=True).items():
            if gid not in native[src]:
                continue
            n_tr += 1
            starts, arr = idx[tgt].get((c, st), ([], []))
            i = bisect.bisect_right(starts, e)
            for j in range(i - 1, max(-1, i - 60), -1):
                bs, be, b = arr[j]
                if be < s:
                    continue
                ov = min(e, be) - max(s, bs) + 1
                if ov <= 0:
                    continue
                if ov >= args.min_frac * (e - s + 1) and ov >= args.min_frac * (be - bs + 1):
                    hit.add((gid, b))
                elif ov >= args.partial_frac * (be - bs + 1) and (be - bs + 1) < args.min_frac * (e - s + 1):
                    partial[gid].add(b)
        print(f"  {src} onto {tgt}: done", file=sys.stderr)
    print(f"{n_tr} primary transfers read", file=sys.stderr)

    # baskets from reciprocal hits
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in genomes:
        for gid in native[g]:
            find(gid)
    recip = 0
    for a, b in hit:
        if (b, a) in hit and a < b:
            recip += 1
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    baskets = defaultdict(list)
    for gid in owner:
        baskets[find(gid)].append(gid)

    # summary
    rows, stats = [], Counter()
    for k, mem in baskets.items():
        per = Counter(owner[m] for m in mem)
        pieces = sorted({p for m in mem for p in partial.get(m, ()) if find(p) != k})
        has_ref = ref in per
        rows.append([k, len(per), max(per.values()), int(has_ref), len(mem), ",".join(sorted(mem)), ",".join(pieces)])
        stats["baskets"] += 1
        stats[f"genomes with a native gene: {'15' if len(per) == len(genomes) else '2-14' if len(per) > 1 else '1'}"] += 1
        stats["with a reference gene"] += has_ref
        stats["with 2+ native genes of one genome"] += max(per.values()) > 1
        stats["with partial members (pieces) attached"] += bool(pieces)
    with open(os.path.join(args.out_dir, "baskets.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["basket", "genomes_with_native", "max_natives_one_genome", "has_reference", "n_natives", "natives",
                    "partial_members"])
        w.writerows(sorted(rows, key=lambda r: (-r[1], r[0])))
    print(f"\nnative genes: {len(owner)}; directed hits: {len(hit)}; reciprocal pairs: {recip}")
    for k in sorted(stats):
        print(f"{stats[k]:7}  {k}")


if __name__ == "__main__":
    main()
