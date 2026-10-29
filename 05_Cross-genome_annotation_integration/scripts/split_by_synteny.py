#!/usr/bin/env python3
"""
Test: do the members of a basket sit at the same place of the genome? A basket is built from sequence (reciprocal
Liftoff transfers); here each member is placed by its neighbours instead. For every member, the baskets of its nearest
--k other-basket genes on each side (in its own genome) are its neighbourhood (orientation ignored). Two members are at
the same place when their neighbourhoods share >= --min-shared baskets. Within a basket, members at the same place are
grouped; the largest group keeps the basket, every other group becomes a basket of its own (a copy elsewhere, e.g. a
paralogue or a translocated gene). Members with fewer than --min-shared neighbours (contig ends) are not judged and stay.
Baskets made only of pieces of another basket are left out (as in choose_models.py).
Writes <out-dir>/baskets_synteny.tsv (same columns as baskets.tsv) and split_members.tsv.
"""
import argparse
import csv
import glob
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_baskets import read_gff   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baskets", required=True)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--min-shared", type=int, default=2)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    rows = list(csv.DictReader(open(a.baskets), delimiter="\t"))
    piece_of = {}
    for r in rows:
        for x in r["partial_members"].split(","):
            if x:
                piece_of.setdefault(x, r["basket"])
    keep = [r for r in rows if not all(m in piece_of and piece_of[m] != r["basket"] for m in r["natives"].split(","))]
    basket_of = {m: r["basket"] for r in keep for m in r["natives"].split(",")}
    owner, pos = {}, {}
    order = defaultdict(list)
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))
        for gid, (c, st, s, e, _) in read_gff(path).items():
            if gid in basket_of:
                owner[gid] = g
                order[(g, c)].append((s, gid))
    for k in order:
        order[k].sort()
        for i, (_, gid) in enumerate(order[k]):
            pos[gid] = (k, i)

    def neighbourhood(gid):
        k, i = pos[gid]
        L, own = order[k], basket_of[gid]
        out = set()
        for step in (-1, 1):
            j, got = i + step, []
            while 0 <= j < len(L) and len(got) < a.k:
                b = basket_of[L[j][1]]
                if b != own and b not in got:
                    got.append(b)
                j += step
            out.update(got)
        return out

    stats = Counter()
    new_rows, split_rows = [], []
    for r in keep:
        mem = [m for m in r["natives"].split(",") if m in pos]
        if len(mem) < 2:
            new_rows.append(r)
            continue
        nb = {m: neighbourhood(m) for m in mem}
        judged = [m for m in mem if len(nb[m]) >= a.min_shared]
        parent = {m: m for m in judged}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i, m1 in enumerate(judged):
            for m2 in judged[i + 1:]:
                if len(nb[m1] & nb[m2]) >= a.min_shared:
                    parent[find(m1)] = find(m2)
        groups = defaultdict(list)
        for m in judged:
            groups[find(m)].append(m)
        groups = sorted(groups.values(), key=lambda g: (-len({owner[x] for x in g}), -len(g)))
        main_g = (groups[0] if groups else []) + [m for m in mem if m not in judged]
        stats["baskets with 2+ members"] += 1
        if len(groups) > 1:
            stats["baskets split"] += 1
            stats["new baskets"] += len(groups) - 1
            for gr in groups[1:]:
                stats["members moved out"] += len(gr)
                stats[f"new basket with members from {'1 genome' if len({owner[x] for x in gr}) == 1 else '2+ genomes'}"] += 1
                nid = sorted(gr)[0] + "_pos"
                d = dict(r)
                d.update(basket=nid, natives=",".join(sorted(gr)), genomes_with_native=len({owner[x] for x in gr}),
                         max_natives_one_genome=max(Counter(owner[x] for x in gr).values()), n_natives=len(gr),
                         partial_members="")
                new_rows.append(d)
                for x in gr:
                    split_rows.append([r["basket"], nid, x, owner[x], pos[x][0][1]])
        d = dict(r)
        d.update(natives=",".join(sorted(main_g)), genomes_with_native=len({owner[x] for x in main_g}),
                 max_natives_one_genome=max(Counter(owner[x] for x in main_g).values()), n_natives=len(main_g))
        new_rows.append(d)
        if Counter(owner[x] for x in main_g).most_common(1)[0][1] > 1:
            stats["kept baskets still with 2+ members of one genome"] += 1
        stats["members not judged (contig end)"] += len(mem) - len(judged)
    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "baskets_synteny.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(new_rows)
    with open(os.path.join(a.out_dir, "split_members.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["basket", "new_basket", "member", "genome", "contig"])
        w.writerows(split_rows)
    print(f"baskets (pieces-only left out): {len(keep)}")
    for k in sorted(stats):
        print(f"  {stats[k]:6}  {k}")


if __name__ == "__main__":
    main()
