#!/usr/bin/env python3
"""Pseudogenes in the pangenes: how many, and are they shared between genomes?

A gene is flagged as pseudogene-like by one of two marks, kept separate:
  biotype     ebi_biotype=pseudogene or gene_biotype=pseudogene (the source annotation calls it a pseudogene;
              Liftoff copies this attribute onto transferred genes)
  broken ORF  valid_ORFs=0 (or valid_ORF=False): a Liftoff transfer with no valid open reading frame in this
              genome. Genes taken straight from a genome's own annotation never went through Liftoff and have
              no such attribute, so they can only be flagged by the biotype mark (the reference genome's own
              genes are in this group).

Gene level: every gene of every genome. Pangene level: only pangenes with exactly one gene per genome, so
that "flagged in a genome" is unambiguous.

Inputs:
  --genes-matrix  pangene_matrix_genes.tr.tab from get_pangenes
  --gff-dir       GFFs whose gene IDs match the matrix, named <genome><--gff-suffix>
Writes <out-prefix>_pangenes.tsv (one row per pangene with at least one flagged gene) and prints a summary.

Standard library only.
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter

from pangenome_pipeline.common.gff import ID_RE


def parse_genes(path):
    """gene ID -> (biotype_pseudo, orf_known, orf_broken)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = f[8]
            i = ID_RE.search(a)
            if not i:
                continue
            biotype = "ebi_biotype=pseudogene" in a or "gene_biotype=pseudogene" in a
            known = broken = False
            m = re.search(r"valid_ORFs=(\d+)", a)
            if m:
                known, broken = True, int(m.group(1)) == 0
            if "valid_ORF=False" in a:
                known, broken = True, True
            elif "valid_ORF=True" in a:
                known = True
            out[i.group(1)] = (biotype, known, broken)
    return out


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "n/a"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_final_annotation_tagged.gff")
    ap.add_argument(
        "--ignore-genomes",
        default="",
        help="comma-separated genomes left out of the pangene level (occupancy is recomputed without them), e.g. the reference, whose own genes cannot be flagged by ORF",
    )
    ap.add_argument("-o", "--out-prefix", required=True)
    args = ap.parse_args()

    with open(args.genes_matrix, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    cols = [(i, g.replace("_for_get_pangenes", "")) for i, g in enumerate(rows[0]) if i > 0 and g.strip()]
    ignore = {g for g in args.ignore_genomes.split(",") if g}
    n_genomes = len([1 for _, g in cols if g not in ignore])

    genes = {}
    for _, genome in cols:
        path = os.path.join(args.gff_dir, genome + args.gff_suffix)
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")
        genes[genome] = parse_genes(path)

    per_genome = {g: Counter() for _, g in cols}
    clean = []
    n_pangenes = n_multi = 0
    for r in rows[1:]:
        if r[0].startswith("chr"):
            continue
        members = []
        for i, genome in cols:
            cell = r[i] if i < len(r) else ""
            if cell not in ("", "-", "NA"):
                members.append((genome, cell.split(",")))
        if not members:
            continue
        n_pangenes += 1
        for genome, ids in members:
            for gid in ids:
                d = genes[genome].get(gid)
                if d is None:
                    continue
                c = per_genome[genome]
                c["genes"] += 1
                c["biotype"] += d[0]
                c["broken"] += d[2]
                c["either"] += d[0] or d[2]
                c["no_orf_info"] += not d[1]
        if any(len(ids) > 1 for _, ids in members):
            n_multi += 1
            continue
        status = {}
        for genome, ids in members:
            if genome in ignore:
                continue
            d = genes[genome].get(ids[0])
            if d is not None:
                status[genome] = d
        if status:
            clean.append((r[0], len(status), status))

    print(
        f"{n_pangenes} pangenes; {len(clean)} with one gene per genome (used for the pangene level), {n_multi} left out\n"
    )
    print("GENE LEVEL - genes in pangenes, per genome")
    print(f"{'genome':11s}{'genes':>7s}{'biotype':>9s}{'broken ORF':>12s}{'either':>8s}{'no ORF info':>13s}")
    tot = Counter()
    for _, g in cols:
        c = per_genome[g]
        tot.update(c)
        print(f"{g:11s}{c['genes']:7d}{c['biotype']:9d}{c['broken']:12d}{c['either']:8d}{c['no_orf_info']:13d}")
    print(
        f"{'ALL':11s}{tot['genes']:7d}{tot['biotype']:9d}{tot['broken']:12d}{tot['either']:8d}{tot['no_orf_info']:13d}"
        f"   ({pct(tot['either'], tot['genes'])} of genes flagged; {pct(tot['no_orf_info'], tot['genes'])} have no ORF information)\n"
    )

    flagged_rows = []
    for pg, occ, status in clean:
        flag = [g for g, d in status.items() if d[0] or d[2]]
        if flag:
            flagged_rows.append(
                (pg, occ, len(flag), sum(status[g][0] for g in flag), sum(status[g][2] for g in flag), ",".join(flag))
            )
    with open(args.out_prefix + "_pangenes.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["pangene", "occupancy", "genomes_flagged", "flagged_by_biotype", "flagged_by_broken_orf", "genomes"]
        )
        w.writerows(flagged_rows)

    n_clean = len(clean)
    n_flag = len(flagged_rows)

    def category(occ):
        if occ >= n_genomes:
            return "Core"
        if occ >= 0.9 * n_genomes:
            return "Soft-core"
        if occ <= 0.15 * n_genomes:
            return "Cloud"
        return "Shell"

    print(
        f"PANGENE LEVEL - {n_clean} pangenes with one gene per genome"
        + (f" (ignoring {', '.join(sorted(ignore))})" if ignore else "")
    )
    print(f"  with at least one flagged gene: {n_flag} ({pct(n_flag, n_clean)})\n")
    print(f"  {'category':10s}{'pangenes':>10s}{'with flagged':>14s}{'%':>8s}")
    for cat in ("Core", "Soft-core", "Shell", "Cloud"):
        tot_c = sum(1 for _, occ, _ in clean if category(occ) == cat)
        fl_c = sum(1 for x in flagged_rows if category(x[1]) == cat)
        if tot_c:
            print(f"  {cat:10s}{tot_c:10d}{fl_c:14d}{pct(fl_c, tot_c):>8s}")

    print("\n  Flagged pangenes: what share of the genomes that have the pangene carry the flagged gene?")
    bins = Counter()
    for x in flagged_rows:
        occ, k = x[1], x[2]
        if occ == 1:
            bins["present in 1 genome only (flagged there = 100%)"] += 1
        elif k == occ:
            bins["all genomes carrying the pangene (2+)"] += 1
        elif k == 1:
            bins["only 1 of the 2+ genomes"] += 1
        elif k / occ >= 0.5:
            bins["2+ genomes, at least half of them"] += 1
        else:
            bins["2+ genomes, less than half of them"] += 1
    for name in (
        "present in 1 genome only (flagged there = 100%)",
        "all genomes carrying the pangene (2+)",
        "2+ genomes, at least half of them",
        "2+ genomes, less than half of them",
        "only 1 of the 2+ genomes",
    ):
        print(f"     {name:52s}{bins[name]:6d}  ({pct(bins[name], n_flag)})")
    multi_occ = [x for x in flagged_rows if x[1] >= 2]
    if multi_occ:
        fr = sorted(x[2] / x[1] for x in multi_occ)
        print(
            f"     among the {len(multi_occ)} flagged pangenes present in 2+ genomes: median share flagged {fr[len(fr) // 2]:.2f}"
        )
    print("\n  Core pangenes (present in every genome considered) by number of genomes flagged:")
    core = Counter(x[2] for x in flagged_rows if x[1] == n_genomes)
    for k in sorted(core):
        print(f"     {k:2d} of {n_genomes}: {core[k]}")
    print(f"\nWrote {args.out_prefix}_pangenes.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
