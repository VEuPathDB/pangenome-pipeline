#!/usr/bin/env python3
"""
The TE line of the explorer: the transposable elements RepeatMasker found in each genome (step 02).

Reads, for every genome, the RepeatMasker table <te-dir>/RM_<genome>/*.out (one folder per genome, named as in
genome_order.txt; e.g. step 02's repeatmasker_contig_filtered/) and writes <out-dir>/te_<genome>.json:
  {"contigs": {contig: [[start, end, strand, name, classification], ...]}}
Simple repeats, low-complexity regions and satellites are left out (--keep-simple keeps them): they are not TEs.

Usage:
    python3 build_te_layer.py --te-dir repeatmasker_contig_filtered --genomes genome_order.txt \\
        -o pangenome_explorer/preliminary/data
"""
import argparse
import glob
import json
import os

SKIP = ("Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA")


def read_out(path, keep_simple):
    by = {}
    for l in open(path):
        f = l.split()
        if len(f) < 11 or not f[0].isdigit():          # the three header lines, blank lines
            continue
        contig, s, e, strand, name, cls = f[4], int(f[5]), int(f[6]), f[8], f[9], f[10]
        if not keep_simple and cls.split("/")[0] in SKIP:
            continue
        by.setdefault(contig, []).append([s, e, "-" if strand == "C" else "+", name, cls])
    for L in by.values():
        L.sort()
    return by


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--te-dir", required=True, help="folder with one RM_<genome>/ folder per genome")
    ap.add_argument("--genomes", required=True, help="genome_order.txt")
    ap.add_argument("--keep-simple", action="store_true")
    ap.add_argument("-o", "--out-dir", required=True)
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    os.makedirs(a.out_dir, exist_ok=True)
    for g in genomes:
        outs = sorted(glob.glob(os.path.join(a.te_dir, f"RM_{g}", "*.out")))
        if not outs:
            print(f"{g}: no RM_{g}/*.out - no TE line")
            continue
        by = read_out(outs[0], a.keep_simple)
        json.dump({"contigs": by}, open(os.path.join(a.out_dir, f"te_{g}.json"), "w"), separators=(",", ":"))
        print(f"{g}: {sum(len(v) for v in by.values())} TEs ({os.path.basename(outs[0])})")


if __name__ == "__main__":
    main()
