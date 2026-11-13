#!/usr/bin/env python3
"""
RNA-seq coverage, from bigWig files, over the CDS of the reference genes added by the pipeline
(not_in_reference_annotation=True in --gff), for the reference gene check (reference_check_job.sh).
Each bigWig file found under --bw-root (searched in every subfolder) that matches --pattern is one column.
Strand: by default ONLY the stranded samples are used, and each gene counts only the reads on its own strand (file
names with "firststrand" = + strand, "secondstrand" = - strand, the same convention as plot_rnaseq_coverage.py), so
the reads of a gene on the other strand (e.g. a neighbour's UTR) do not count. --all-samples adds the unstranded
samples (both strands). When a sample folder has strand files and also the combined file of the same reads, the
combined file is never used.

Per item: all_samples = the mean coverage summed over all samples; the best sample, its mean and the fraction of bases
with reads there; how many samples have reads; the profile = upstream flank, 10 equal parts, downstream flank (mean
of the samples). For comparison, --baseline reference genes of the reference's own annotation (genes with an
original_ID, picked at random, always the same ones) are measured the same way: percentile_own = the share of them
with less coverage than the item.

Output (--out, TSV): a first line "# baseline<TAB>median_mean_coverage<TAB><value>" (the median, over the own genes
measured, of their mean coverage per sample: the level of a typical reference gene, in the profile's unit), then kind, id, contig, start, end, all_samples, percentile_own, best_sample, best_mean, best_breadth,
samples_with_reads, profile (12 values, comma-separated)

Needs pyBigWig.
"""
import argparse
import bisect
import fnmatch
import gzip
import os
import random
import re
import sys

import pyBigWig

att = lambda s, k: (re.search(r"(?:^|;)" + k + r"=([^;\n]*)", s) or [None, ""])[1]
openf = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_gff(gff):
    """gene -> (contig, start, end, strand, new, own); gene -> CDS segments (0-based, half-open)"""
    genes, t2g, cds = {}, {}, {}
    for l in openf(gff):
        if l.startswith("#") or l.count("\t") < 8:
            continue
        f = l.rstrip("\n").split("\t")
        if f[2] == "gene":
            genes[att(f[8], "ID")] = (f[0], int(f[3]), int(f[4]), f[6],
                                      att(f[8], "not_in_reference_annotation") == "True", bool(att(f[8], "original_ID")))
        elif f[2] in ("mRNA", "transcript"):
            t2g[att(f[8], "ID")] = att(f[8], "Parent")
        elif f[2] == "CDS":
            for p in att(f[8], "Parent").split(","):
                cds.setdefault(p, []).append((f[0], int(f[3]) - 1, int(f[4])))
    per = {}
    for t, segs in cds.items():
        per.setdefault(t2g.get(t, t), set()).update(segs)
    return genes, {g: sorted(s, key=lambda x: x[1]) for g, s in per.items()}


def cov(bw, chroms, segs):
    """(mean, fraction of bases with reads) over the segments; a contig missing from the bigWig = no reads"""
    tot = covd = n = 0
    for c, s, e in segs:
        s = max(0, s)
        if c in chroms:
            e2 = min(e, chroms[c])
            for a, b, v in (bw.intervals(c, s, e2) if e2 > s else None) or []:
                a, b = max(a, s), min(b, e2)
                if b > a and v > 0:
                    tot += (b - a) * v
                    covd += b - a
        n += max(0, e - s)
    return (tot / n, covd / n) if n else (0.0, 0.0)


def parts(segs, strand, flank, nb=10):
    """upstream flank, nb equal parts of the bases (5' to 3'), downstream flank"""
    c = segs[0][0]
    bases = sum(e - s for _, s, e in segs)
    cuts = [round(i * bases / nb) for i in range(nb + 1)]
    order = segs if strand != "-" else segs[::-1]
    out, pos = [[] for _ in range(nb)], 0
    for _, s, e in order:
        n = e - s
        for i in range(nb):
            a, b = max(cuts[i], pos), min(cuts[i + 1], pos + n)
            if b > a:
                out[i].append((c, s + a - pos, s + b - pos) if strand != "-" else (c, e - (b - pos), e - (a - pos)))
        pos += n
    lo, hi = segs[0][1], segs[-1][2]
    left, right = [(c, max(0, lo - flank), lo)], [(c, hi, hi + flank)]
    return ([right] if strand == "-" else [left]) + out + ([left] if strand == "-" else [right])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, help="the reference's final GFF (step 05)")
    ap.add_argument("--bw-root", required=True, help="folder searched (with its subfolders) for bigWig files")
    ap.add_argument("--pattern", default="unique_*unlogged*.bw", help="bigWig file names to use")
    ap.add_argument("--baseline", type=int, default=500, help="own reference genes measured for the percentile")
    ap.add_argument("--flank", type=int, default=300)
    ap.add_argument("--all-samples", action="store_true", help="also the unstranded samples (reads of both strands)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes, segs = read_gff(a.gff)
    items = [("new_gene", g, x[0], x[1], x[2], x[3], segs.get(g) or [(x[0], x[1] - 1, x[2])])
             for g, x in genes.items() if x[4]]
    own = sorted(g for g, x in genes.items() if x[5] and segs.get(g))
    base = random.Random(1).sample(own, min(a.baseline, len(own)))
    bws = []
    for d, _, fs in os.walk(a.bw_root):
        fs = [f for f in fs if fnmatch.fnmatch(f, a.pattern)]
        stranded = [f for f in fs if "firststrand" in f or "secondstrand" in f]
        bws += [os.path.join(d, f) for f in (stranded or ([] if not a.all_samples else fs))]   # combined file: only without strand files
    bws.sort()
    fstrand = lambda p: "+" if "firststrand" in os.path.basename(p) else "-" if "secondstrand" in os.path.basename(p) else None
    print(f"{len(items)} new genes, {len(base)} own genes for the percentile, {len(bws)} bigWig files "
          f"({sum(1 for p in bws if fstrand(p) == '+')} + strand, {sum(1 for p in bws if fstrand(p) == '-')} - strand, "
          f"{sum(1 for p in bws if fstrand(p) is None)} unstranded); each gene counts only its own strand", file=sys.stderr)
    if not bws:
        sys.exit(f"ERROR: no bigWig file matching {a.pattern} under {a.bw_root}")

    names = [os.path.relpath(p, a.bw_root) for p in bws]
    res = [dict(sum=0.0, best=("", 0.0, 0.0), reads=0, prof=[0.0] * 12, n=0) for _ in items]
    pp = [parts(it[6], it[5], a.flank) for it in items]
    bsum, bn = [0.0] * len(base), [0] * len(base)
    for j, p in enumerate(bws):
        bw = pyBigWig.open(p)
        ch = bw.chroms()
        fs_ = fstrand(p)
        for k, it in enumerate(items):
            if fs_ is not None and fs_ != it[5]:
                continue                           # a stranded file of the other strand
            res[k]["n"] += 1
            m, br = cov(bw, ch, it[6])
            r = res[k]
            r["sum"] += m
            r["reads"] += m > 0
            if m > r["best"][1]:
                r["best"] = (names[j], m, br)
            for i, ps in enumerate(pp[k]):
                r["prof"][i] += cov(bw, ch, ps)[0]
        for k, g in enumerate(base):
            if fs_ is not None and fs_ != genes[g][3]:
                continue
            bsum[k] += cov(bw, ch, segs[g])[0]
            bn[k] += 1
        bw.close()
        print(f"[{j + 1}/{len(bws)}] {names[j]}", file=sys.stderr)

    per = sorted(b / n for b, n in zip(bsum, bn) if n)   # mean per sample of each own gene
    bsum.sort()
    with open(a.out, "w") as fh:
        med = per[len(per) // 2] if per else 0.0
        fh.write(f"# baseline\tmedian_mean_coverage\t{med:.3f}\n")
        fh.write("kind\tid\tcontig\tstart\tend\tall_samples\tpercentile_own\tbest_sample\tbest_mean\tbest_breadth\t"
                 "samples_with_reads\tprofile\n")
        for it, r in zip(items, res):
            pct = f"{100 * bisect.bisect_left(bsum, r['sum']) / len(bsum):.1f}" if bsum else ""
            fh.write("\t".join([it[0], it[1], it[2], str(it[3]), str(it[4]), f"{r['sum']:.2f}", pct, r["best"][0],
                                f"{r['best'][1]:.2f}", f"{r['best'][2]:.3f}", str(r["reads"]),
                                ",".join(f"{v / max(r['n'], 1):.2f}" for v in r["prof"])]) + "\n")
    print("->", a.out, file=sys.stderr)


if __name__ == "__main__":
    main()
