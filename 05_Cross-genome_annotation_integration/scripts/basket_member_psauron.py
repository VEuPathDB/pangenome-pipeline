#!/usr/bin/env python3
"""
PSAURON of every member of selected baskets, to evaluate the winner rule (nothing is changed).

Selected baskets:
  multi_winner       the basket has two or more winners (often a reference gene plus another strain's model next to
                     it, while the strains' own genes cover both)
  ref_invalid        a reference winner whose model has a valid ORF in at most --max-valid-share of the genomes where
                     it was placed (its transfers are broken in most strains)
A basket can be both.

For each member (every native gene of the basket, in every genome): its genome, whether it won, whether it is a
reference gene, its CDS length (first transcript of the native GFF), whether its own ORF is invalid (step 05's
decide/own_invalid_orf.tsv) and its PSAURON in-frame score (step 04, keyed by gene or transcript ID).

Outputs (in --out-dir):
  basket_members_psauron.tsv   one row per member
  basket_psauron_summary.tsv   one row per basket: the winners and their scores, the best-scoring member that did not
                               win, and better_member=yes when that member scores more than --margin above the best
                               winner (or the winners have no score) and is at least as long as the longest winner
  summary.tsv                  counts

Usage:
    python3 basket_member_psauron.py --basket-dir basket_annotation --native-dir basket_annotation/input/native \\
        --psauron-dir psauron --out-dir basket_annotation/psauron_members
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

from baskets_v2 import read_psauron

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def read_natives(native_dir):
    """gene ID -> (genome, CDS length of its first transcript); transcript ID -> gene ID."""
    genes, tx2gene = {}, {}
    paths = sorted(glob.glob(os.path.join(native_dir, "*.gff3")) + glob.glob(os.path.join(native_dir, "*.gff")))
    if not paths:
        sys.exit(f"ERROR: no .gff3/.gff in {native_dir}")
    for path in paths:
        genome = re.sub(r"\.gff3?$", "", os.path.basename(path))
        ids, first, cds = [], {}, Counter()
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
            if i and not p:                                 # a gene, whatever its type (gene, protein_coding_gene...)
                ids.append(i.group(1))
            elif f[2] in ("mRNA", "transcript", "pseudogenic_transcript") and i and p:
                tx2gene[i.group(1)] = p.group(1)
                first.setdefault(p.group(1), i.group(1))
            elif f[2] == "CDS" and p:
                cds[p.group(1).split(",")[0]] += int(f[4]) - int(f[3]) + 1
        for g in ids:
            genes[g] = (genome, cds.get(first.get(g), 0))
    return genes, tx2gene


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--basket-dir", required=True, help="step 05 output (baskets_v2.tsv, winners.tsv, decide/)")
    ap.add_argument("--native-dir", required=True, help="the native GFFs step 05 read (basket_annotation/input/native)")
    ap.add_argument("--psauron-dir", required=True, help="step 04 PSAURON scores (*.csv)")
    ap.add_argument("--max-valid-share", type=float, default=0.5,
                    help="ref_invalid: reference winner valid in at most this share of the genomes it is placed in")
    ap.add_argument("--margin", type=float, default=0.05, help="better_member: score above the best winner by this")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    genes, tx2gene = read_natives(a.native_dir)
    score = {}
    for f in sorted(glob.glob(os.path.join(a.psauron_dir, "*.csv"))):
        for k, v in read_psauron(f).items():
            gk = k if k in genes else tx2gene.get(k, tx2gene.get(re.sub(r"-p\d+$", "", k), k))
            score[gk] = v
    print(f"{len(genes)} native genes, {len(score)} PSAURON scores", file=sys.stderr)
    win = {r["gene"]: r for r in csv.DictReader(open(os.path.join(a.basket_dir, "winners.tsv")), delimiter="\t")}
    own_bad = set()
    p = os.path.join(a.basket_dir, "decide", "own_invalid_orf.tsv")
    if os.path.exists(p):
        own_bad = {r["gene"] for r in csv.DictReader(open(p), delimiter="\t")}

    def ref_invalid(w):
        r = win.get(w, {})
        try:
            return r["is_reference"] == "1" and int(r["valid_orf_genomes"]) <= a.max_valid_share * int(r["present_in"])
        except (KeyError, ValueError):
            return False

    members_out, summ = [], []
    for r in csv.DictReader(open(os.path.join(a.basket_dir, "baskets_v2.tsv")), delimiter="\t"):
        ws = [w for w in r["winners"].split(";") if w]
        cat = (["multi_winner"] if len(ws) >= 2 else []) + (["ref_invalid"] if any(ref_invalid(w) for w in ws) else [])
        if not cat:
            continue
        cat = "+".join(cat)
        ms = [m for m in r["members"].split(",") if m]
        for m in ms:
            g, cds = genes.get(m, ("?", 0))
            sc = score.get(m)
            members_out.append([r["winners"], cat, m, g, "yes" if m in ws else "no",
                                "yes" if win.get(m, {}).get("is_reference") == "1" else "no", cds,
                                "invalid" if m in own_bad else "", f"{sc[1]:.3f}" if sc else "", sc[0] if sc else ""])
        wsc = [score[w][1] for w in ws if w in score]
        others = [m for m in ms if m not in ws and m in score]
        best = max(others, key=lambda m: score[m][1]) if others else ""
        bsc = score[best][1] if best else None
        wlen = max([genes.get(w, ("", 0))[1] for w in ws] or [0])
        better = bool(best) and (not wsc or bsc > max(wsc) + a.margin) and genes.get(best, ("", 0))[1] >= wlen
        summ.append([r["winners"], cat, len(ms), sum(1 for m in ms if m in score),
                     sum(1 for m in ms if m in score and score[m][1] >= 0.9),
                     ";".join(f"{w}:{score[w][1]:.3f}" if w in score else f"{w}:-" for w in ws),
                     ";".join(f"{w}:{win.get(w, {}).get('valid_orf_genomes', '')}/{win.get(w, {}).get('present_in', '')}"
                              for w in ws),
                     best, genes.get(best, ("", ""))[0] if best else "", f"{bsc:.3f}" if best else "",
                     genes.get(best, ("", ""))[1] if best else "", wlen, "yes" if better else "no"])

    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "basket_members_psauron.tsv"), "w") as fh:
        fh.write("winners\tcategory\tmember\tgenome\tis_winner\tis_reference\tcds_bp\town_orf\tpsauron\t"
                 "psauron_is_protein\n")
        for row in members_out:
            fh.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(a.out_dir, "basket_psauron_summary.tsv"), "w") as fh:
        fh.write("winners\tcategory\tmembers\tmembers_scored\tmembers_psauron_0.9_or_more\twinner_psauron\t"
                 "winner_valid_orf(genomes/placed)\tbest_other_member\tits_genome\tits_psauron\tits_cds_bp\t"
                 "longest_winner_cds_bp\tbetter_member\n")
        for row in sorted(summ, key=lambda x: (x[12] != "yes", x[1], x[0])):
            fh.write("\t".join(map(str, row)) + "\n")
    c = Counter(x[1] for x in summ)
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        fh.write(f"baskets\t{len(summ)}\n")
        for k, n in sorted(c.items()):
            fh.write(f"baskets_{k}\t{n}\n")
        fh.write(f"baskets_with_a_better_member\t{sum(1 for x in summ if x[12] == 'yes')}\n")
        fh.write(f"members\t{len(members_out)}\nmembers_without_psauron\t{sum(1 for x in members_out if not x[8])}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read())


if __name__ == "__main__":
    main()
