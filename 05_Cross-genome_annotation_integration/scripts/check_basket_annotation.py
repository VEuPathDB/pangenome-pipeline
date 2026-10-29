#!/usr/bin/env python3
"""
Checks of the GFFs written by build_basket_annotation.py:
  - every gene ID is unique across all genomes
  - every gene has a CDS (except pseudogenes), and its transcript, exons and CDS lie inside it, on its contig and strand
  - no two genes of a genome lie on the same strand with CDS spans overlapping >= 50% of the smaller (stacked genes)
  - every protein-coding gene of the reference annotation is written in the reference genome, unless it is listed in
    removed_panel_genes.tsv
Optionally (--cases), prints what each listed genome got for each listed panel gene (basket).
Usage: python3 check_basket_annotation.py --gff-dir gff_similar --genomes genome_order.txt --reference-gff ref.gff3
           [--cases "GENOME:PANEL_GENE,GENOME:PANEL_GENE"]
"""
import argparse
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_baskets import read_gff   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--reference-gff", required=True, help="the reference genome's own annotation")
    ap.add_argument("--suffix", default="_basket_annotation.gff")
    ap.add_argument("--cases", default="")
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    ref = genomes[0]
    ref_genes = set(read_gff(a.reference_gff))
    removed = set()
    rp = os.path.join(a.gff_dir, "removed_panel_genes.tsv")
    if os.path.exists(rp):
        removed = {l.split("\t")[0] for l in open(rp).read().splitlines()[1:]}
    ids, problems, seen_ref, per = Counter(), Counter(), set(), {}
    for g in genomes:
        genes, kids, cds = {}, defaultdict(list), defaultdict(list)
        for line in open(os.path.join(a.gff_dir, g + a.suffix)):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            i = re.search(r"ID=([^;]+)", f[8]).group(1)
            ids[i] += 1
            if f[2] == "gene":
                genes[i] = f
                per[(g, re.search(r"basket=([^;]+)", f[8]).group(1))] = f[8]
                if g == ref:
                    seen_ref.add(re.search(r"model_from=([^;]+)", f[8]).group(1))
            else:
                gid = re.search(r"Parent=([^;]+)", f[8]).group(1).split(".t")[0]
                kids[gid].append(f)
                if f[2] == "CDS":
                    cds[gid].append((int(f[3]), int(f[4])))
        for i, f in genes.items():
            if not cds.get(i) and "pseudogene=True" not in f[8]:   # pseudogenes: gene + transcript + exon
                problems[f"{g}: gene without CDS"] += 1
            if any(k[0] != f[0] or k[6] != f[6] or int(k[3]) < int(f[3]) or int(k[4]) > int(f[4]) for k in kids[i]):
                problems[f"{g}: feature outside its gene"] += 1
        L = sorted((genes[i][0], genes[i][6], min(x for x, _ in c), max(y for _, y in c)) for i, c in cds.items())
        for x in range(len(L)):
            for y in range(x + 1, min(x + 8, len(L))):
                c1, t1, s1, e1 = L[x]
                c2, t2, s2, e2 = L[y]
                if c1 == c2 and t1 == t2 and min(e1, e2) - max(s1, s2) + 1 >= 0.5 * min(e1 - s1 + 1, e2 - s2 + 1):
                    problems[f"{g}: stacked genes (same strand, >= 50%)"] += 1
    dup = sum(1 for v in ids.values() if v > 1)
    missing = sorted(ref_genes - seen_ref - removed)
    print(f"duplicated IDs: {dup}")
    print(f"reference genes written in {ref}: {len(ref_genes & seen_ref)} of {len(ref_genes)}; "
          f"removed on purpose: {len(ref_genes & removed)}; missing: {len(missing)} {missing[:10]}")
    for k, v in sorted(problems.items()):
        print(f"PROBLEM {v:5}  {k}")
    if not problems and not dup and not missing:
        print("all checks passed")
    for case in [c for c in a.cases.split(",") if c]:
        g, b = case.split(":")
        attrs = per.get((g, b))
        if attrs:
            keep = re.findall(r"(model_origin=[^;]+|model_from=[^;]+|orf_invalid=True|low_coverage_transfer=True|selenoprotein=True)", attrs)
            print(f"  {g:10} {b:18} {' '.join(keep)}")
        else:
            print(f"  {g:10} {b:18} (no gene)")


if __name__ == "__main__":
    main()
