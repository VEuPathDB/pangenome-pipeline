#!/usr/bin/env python3
"""
Apply only the safe fixes to the get_pangenes matrix, and explain every absence that is left.

Safe fix: a gene that audit_absences.py found at the expected place of a pangene P, with the same sequence as P's genes,
but inside another pangene Q that has no genome in common with P (code 1a, same_gene_movable). get_pangenes split one
gene into two pangenes; the gene is moved from Q to P. Only the matrix changes: no gene is added to or removed from the
GFFs. A move is skipped (left for inspection) when:
  - the gene is claimed by more than one pangene;
  - P and Q claim genes of each other: the genes go to the pangene with more genomes, the other direction is skipped;
  - the pangene is both a source and a target of moves (chains).
A pangene left with no gene is removed.

Safe join (after the moves): step 05's basket= tags say which genes are the same gene. When get_pangenes split the genes
of one basket over several pangenes, and every one of these pangenes holds only genes of that basket and no genome has
two of them, the pangenes are joined into the one with the most genomes. Splits where a pangene also holds genes of
other baskets, or a genome has two genes of the basket, are not joined (basket_splits_not_joined.tsv).

Every absence left (pangene P with no gene in genome g, after the moves) gets a verdict and a reason in plain words:
  verdict  absent          the gene is not in that strain (2a other gene in place, 4 place shorter, 5 empty place,
                           8 region absent)
           assembly        cannot tell: the place is at a contig edge or the neighbours are split (6, 7)
           copy_dropped    a Liftoff transfer of P is there, but step 05 did not use it (2b, 3); column step05 says why
           elsewhere       the gene is there but in another pangene, not moved (1a skipped, 1b, 1c, gene moved out)
and the column step05 says what step 05 did with P's panel genes (basket= tags) in that genome: placed elsewhere, dropped
for another gene at the same place (stacked_dropped.tsv), or left out (absent.tsv, with its reason).

Inputs: the get_pangenes matrix (pangene_matrix_genes.tr.tab, pangenes in rows; or pangene_matrix_genes.tab, genomes in
rows), audit_absences.py's absences.tsv, the step 05 final GFFs (for basket=), absent.tsv and stacked_dropped.tsv.
Outputs (in --out-dir):
  pangene_matrix_genes.tr.tab  the fixed matrix (pangenes in rows), and pangene_matrix_genes.tab (genomes in rows)
  moves.tsv                    one row per gene moved: genome, gene, from, to, identity, coverage
  skipped_moves.tsv            1a cases not moved, and why
  joins.tsv                    one row per pangene joined into another: basket, kept pangene, joined pangene, genomes
  basket_splits_not_joined.tsv baskets still split over several pangenes, and why they were not joined
  absences_explained.tsv       one row per absence left: verdict, reason, step 05's reason, and the audit's details
  summary.tsv                  counts

Usage:
    python3 fix_pangenome_safe.py --genes-matrix <pangenes>/pangene_matrix_genes.tr.tab \\
        --absences pangenome_inspection/absence_audit/absences.tsv \\
        --gff-dir basket_annotation/gff_prefer_panel_model_final \\
        --absent basket_annotation/decide/absent.tsv \\
        --stacked-dropped basket_annotation/gff_prefer_panel_model/stacked_dropped.tsv \\
        --out-dir pangenome_fixed
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

EMPTY = {"", "-", "NA"}
VERDICT = {"1a": "elsewhere", "1b": "elsewhere", "1c": "elsewhere", "2a": "absent", "2b": "copy_dropped",
           "3": "copy_dropped", "4": "absent", "5": "absent", "6": "assembly", "7": "assembly", "8": "absent"}
REASON = {
    "1a": "the gene is here, in another pangene (match_pangene); not moved, see skipped_moves.tsv",
    "1b": "the gene is here, in another pangene (match_pangene) that also has genes in strains with this pangene "
          "(paralog or duplication?)",
    "1c": "the gene is here, but in no pangene (left out by get_pangenes, e.g. CDS < 100 bp)",
    "2a": "another gene sits in the place, and there is no Liftoff transfer of this gene",
    "2b": "another gene sits in the place, over a Liftoff transfer of this gene (step 05 kept the other gene)",
    "3": "nothing in the place but a Liftoff transfer of this gene, not used by step 05",
    "4": "nothing in the place, and the place is much shorter than in the other strains (deletion?)",
    "5": "nothing in the place, no transfer: the gene is missing",
    "6": "the place is at a contig edge: cannot tell",
    "7": "the neighbours are on different contigs or far apart (assembly break or rearrangement): cannot tell",
    "8": "the region is not found in this strain",
}


def read_matrix(path):
    """-> (header, genomes, [pangene], {pangene: {genome: [genes]}}, extra rows kept as they are)."""
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    if not path.endswith(".tr.tab"):                     # genomes in rows: transpose
        n = max(len(r) for r in rows)
        rows = [list(c) for c in zip(*[r + [""] * (n - len(r)) for r in rows])]
    header = rows[0]
    genomes = header[1:]
    order, cells, extra = [], {}, []
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        if r[0].startswith("chr:"):
            extra.append(r)
            continue
        order.append(r[0])
        cells[r[0]] = {g: [x for x in re.split(r"[,;]", c) if x] for g, c in zip(genomes, r[1:]) if c not in EMPTY}
    return header, genomes, order, cells, extra


def write_matrix(path, header, genomes, order, cells, extra):
    rows = [header] + extra
    for p in order:
        rows.append([p] + [",".join(cells[p][g]) if g in cells[p] else "-" if g else "" for g in genomes])
    with open(path + ".tr.tab", "w") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    n = max(len(r) for r in rows)
    with open(path + ".tab", "w") as fh:
        for c in zip(*[r + [""] * (n - len(r)) for r in rows]):
            fh.write("\t".join(c) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--absences", required=True, help="audit_absences.py's absences.tsv")
    ap.add_argument("--gff-dir", required=True, help="step 05 final GFFs (basket= tags)")
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--absent", required=True, help="step 05's decide/absent.tsv")
    ap.add_argument("--stacked-dropped", required=True, help="step 05's stacked_dropped.tsv")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    header, genomes, order, cells, extra = read_matrix(a.genes_matrix)
    col = {g.replace("_for_get_pangenes", ""): g for g in genomes}    # absences.tsv uses the short genome names
    print(f"{len(order)} pangenes, {len(genomes)} genomes", file=sys.stderr)

    # 1. the safe moves
    absences = list(csv.DictReader(open(a.absences), delimiter="\t"))
    cand = [r for r in absences if r["code"] == "1a"]
    size = {p: len(c) for p, c in cells.items()}
    claims = Counter(r["match"] for r in cand)
    pairs = {(r["pangene"], r["match_pangene"]) for r in cand}
    moves, skipped = [], []
    for r in cand:
        P, Q, x = r["pangene"], r["match_pangene"], r["match"]
        why = ""
        if claims[x] > 1:
            why = "gene claimed by more than one pangene"
        elif (Q, P) in pairs and (size.get(Q, 0), Q) > (size.get(P, 0), P):
            why = f"{P} and {Q} claim genes of each other: genes go to {Q}, the pangene with more genomes"
        elif P not in cells or Q not in cells:
            why = "pangene not in the matrix"
        if why:
            skipped.append((r, why))
        else:
            moves.append(r)
    # chains: a pangene that gives genes and also receives genes
    targets, sources = {r["pangene"] for r in moves}, {r["match_pangene"] for r in moves}
    chain = targets & sources
    for r in [m for m in moves if m["pangene"] in chain or m["match_pangene"] in chain]:
        moves.remove(r)
        skipped.append((r, "pangene both gives and receives genes (chain)"))

    moved_out = defaultdict(dict)          # (Q, genome) -> P
    done = []
    for r in moves:
        P, Q, x, g = r["pangene"], r["match_pangene"], r["match"], col.get(r["genome"], r["genome"])
        if g in cells[P] or x not in cells[Q].get(g, []):
            skipped.append((r, "matrix cell changed by another move"))
            continue
        cells[Q][g].remove(x)
        if not cells[Q][g]:
            del cells[Q][g]
            moved_out[Q][g] = P
        cells[P][g] = [x]
        done.append(r)
    removed = [p for p in order if not cells[p]]
    order = [p for p in order if cells[p]]

    basket, where = {}, defaultdict(list)
    for path in glob.glob(os.path.join(a.gff_dir, "*" + a.gff_suffix)):
        g = os.path.basename(path)[:-len(a.gff_suffix)]
        for line in open(path):
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                gid = re.search(r"(?:^|;)ID=([^;]+)", f[8]).group(1)
                b = re.search(r"basket=([^;\n]+)", f[8])
                if b:
                    basket[gid] = b.group(1)
                    where[(b.group(1), g)].append(gid)

    # 1b. the safe joins: one basket split over pangenes that hold only that basket, one gene per genome
    of_basket = defaultdict(set)                        # basket -> the pangenes with a gene of it
    for p in order:
        for xs in cells[p].values():
            for x in xs:
                of_basket[basket.get(x)].add(p)
    holds = {p: {basket.get(x) for xs in cells[p].values() for x in xs} for p in order}
    joins, not_joined, joined_away = [], [], set()
    for b, ps in sorted((b, ps) for b, ps in of_basket.items() if b and len(ps) > 1):
        ps = sorted(ps)
        n_genomes = Counter(g for p in ps for g, xs in cells[p].items() for x in xs if basket.get(x) == b)
        if any(holds[p] != {b} for p in ps):
            not_joined.append((b, ps, "a pangene also holds genes of other baskets: " + ",".join(
                sorted({o or "no basket" for p in ps for o in holds[p]} - {b}))))
        elif max(n_genomes.values()) > 1:
            not_joined.append((b, ps, "a genome has two genes of the basket: " + ",".join(
                sorted(g.replace("_for_get_pangenes", "") for g, n in n_genomes.items() if n > 1))))
        else:
            keep = max(ps, key=lambda p: (len(cells[p]), p))
            for p in ps:
                if p == keep:
                    continue
                joins.append((b, keep, p, sorted(cells[p])))
                cells[keep].update(cells[p])
                cells[p] = {}
                joined_away.add(p)
    order = [p for p in order if cells[p]]
    filled = {(keep, g.replace("_for_get_pangenes", "")) for _, keep, _, gs in joins for g in gs}

    # 2. what step 05 did with each pangene's panel genes
    pang_of = {x: p for p in order for xs in cells[p].values() for x in xs}
    absent = {(r["panel_gene"], r["genome"]): r["reason"] for r in csv.DictReader(open(a.absent), delimiter="\t")}
    stacked = {(r["panel_gene_dropped"], r["genome"]): r["kept_panel_gene"]
               for r in csv.DictReader(open(a.stacked_dropped), delimiter="\t")}

    def step05(P, g):
        out = []
        for b in sorted({basket.get(x) for xs in cells.get(P, {}).values() for x in xs} - {None}):
            if where.get((b, g)):
                out.append(f"{b}: placed as " + ",".join(f"{x} (in {pang_of.get(x, 'no pangene')})"
                                                        for x in where[(b, g)]))
            elif (b, g) in stacked:
                out.append(f"{b}: dropped, the place was kept for {stacked[(b, g)]}")
            elif (b, g) in absent:
                out.append(f"{b}: left out ({absent[(b, g)]})")
            else:
                out.append(f"{b}: no record")
        return "; ".join(out)

    # 3. the absences left
    fixed = {(r["pangene"], r["genome"]) for r in done}
    cols = [c for c in absences[0].keys() if c not in ("pangene", "genome", "code", "category")]
    out = []
    for r in absences:
        if (r["pangene"], r["genome"]) in fixed or (r["pangene"], r["genome"]) in filled \
                or r["pangene"] in removed or r["pangene"] in joined_away:
            continue
        out.append([r["pangene"], r["genome"], VERDICT[r["code"]], r["code"], REASON[r["code"]],
                    step05(r["pangene"], r["genome"])] + [r[c] for c in cols])
    for Q, gs in moved_out.items():
        if Q in removed or Q in joined_away:
            continue
        for g, P in gs.items():
            gn = g.replace("_for_get_pangenes", "")
            out.append([Q, gn, "elsewhere", "moved", f"the gene was moved to {P} (same gene, split by get_pangenes)",
                        step05(Q, gn)] + [""] * len(cols))

    os.makedirs(a.out_dir, exist_ok=True)
    write_matrix(os.path.join(a.out_dir, "pangene_matrix_genes"), header, genomes, order, cells, extra)
    with open(os.path.join(a.out_dir, "moves.tsv"), "w") as fh:
        fh.write("genome\tgene\tfrom_pangene\tto_pangene\tidentity\tcoverage\n")
        for r in done:
            fh.write(f"{r['genome']}\t{r['match']}\t{r['match_pangene']}\t{r['pangene']}\t{r['match_identity']}\t"
                     f"{r['match_cov']}\n")
    with open(os.path.join(a.out_dir, "skipped_moves.tsv"), "w") as fh:
        fh.write("genome\tgene\tin_pangene\tclaimed_by_pangene\twhy_not_moved\n")
        for r, why in skipped:
            fh.write(f"{r['genome']}\t{r['match']}\t{r['match_pangene']}\t{r['pangene']}\t{why}\n")
    with open(os.path.join(a.out_dir, "joins.tsv"), "w") as fh:
        fh.write("basket\tkept_pangene\tjoined_pangene\tgenomes\n")
        for b, keep, p, gs in joins:
            fh.write(f"{b}\t{keep}\t{p}\t{','.join(g.replace('_for_get_pangenes', '') for g in gs)}\n")
    with open(os.path.join(a.out_dir, "basket_splits_not_joined.tsv"), "w") as fh:
        fh.write("basket\tpangenes\twhy_not_joined\n")
        for b, ps, why in not_joined:
            fh.write(f"{b}\t{','.join(ps)}\t{why}\n")
    with open(os.path.join(a.out_dir, "absences_explained.tsv"), "w") as fh:
        fh.write("\t".join(["pangene", "genome", "verdict", "audit_code", "reason", "step05"] + cols) + "\n")
        for r in sorted(out):
            fh.write("\t".join(map(str, r)) + "\n")
    vc = Counter(r[2] for r in out)
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        fh.write(f"pangenes_before\t{len(order) + len(removed) + len(joined_away)}\npangenes_after\t{len(order)}\n"
                 f"genes_moved\t{len(done)}\nmoves_skipped\t{len(skipped)}\npangenes_removed\t{len(removed)}\n"
                 f"baskets_joined\t{len({j[0] for j in joins})}\npangenes_joined_away\t{len(joined_away)}\n"
                 f"basket_splits_not_joined\t{len(not_joined)}\n"
                 f"absences_before\t{len(absences)}\nabsences_after\t{len(out)}\n")
        for v in ("absent", "assembly", "copy_dropped", "elsewhere"):
            fh.write(f"verdict_{v}\t{vc[v]}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read())


if __name__ == "__main__":
    main()
