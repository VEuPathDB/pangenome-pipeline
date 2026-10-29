#!/usr/bin/env python3
"""
Remove whole panel genes (e.g. transposable elements) from the final GFFs and from the pangene matrix.

A panel gene is a step 05 basket: every gene tagged basket=<panel gene> in any genome is removed, with all its features
(mRNA/transcript, exon, CDS...: every line whose Parent chain leads to the gene). The other lines of the GFFs are written
unchanged. In the matrix the removed genes leave their pangene; a pangene left with no gene is removed.

--remove is a TSV with the columns panel_gene and reason (one panel gene per row, '#' lines ignored).

Outputs (in --out-dir):
  gff/<genome><suffix>             the GFFs without the removed genes
  pangene_matrix_genes.tr.tab/.tab the matrix without them (pangenes in rows / genomes in rows)
  removed_genes.tsv                one row per gene removed: genome, gene, panel gene, location, pangene, reason
  absences_explained.tsv           (with --absences) the same table without the rows of the pangenes removed
  summary.tsv                      counts

Usage:
    python3 remove_panel_genes.py --remove genes_to_remove.tsv \\
        --gff-dir basket_annotation/gff_prefer_panel_model_final \\
        --genes-matrix pangenome_fixed/pangene_matrix_genes.tr.tab \\
        --absences pangenome_fixed/absences_explained.tsv --out-dir pangenome_final
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter

from fix_pangenome_safe import read_matrix, write_matrix

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
BASKET_RE = re.compile(r"(?:^|;)basket=([^;\n]+)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--remove", required=True, help="TSV: panel_gene, reason")
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--absences", default="", help="absences_explained.tsv of fix_pangenome_safe.py (optional)")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    reason = {r["panel_gene"]: r["reason"] for r in csv.DictReader(
        (l for l in open(a.remove) if not l.startswith("#")), delimiter="\t")}
    print(f"{len(reason)} panel genes to remove", file=sys.stderr)

    # 1. GFFs
    os.makedirs(os.path.join(a.out_dir, "gff"), exist_ok=True)
    removed = []                                   # (genome, gene, panel gene, location)
    for path in sorted(glob.glob(os.path.join(a.gff_dir, "*" + a.gff_suffix))):
        genome = os.path.basename(path)[:-len(a.gff_suffix)]
        lines = open(path).readlines()
        drop = set()
        for line in lines:                         # the genes of the panel genes to remove
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene":
                b = BASKET_RE.search(f[8])
                if b and b.group(1) in reason:
                    gid = ID_RE.search(f[8]).group(1)
                    drop.add(gid)
                    removed.append((genome, gid, b.group(1), f"{f[0]}:{f[3]}-{f[4]}{f[6]}"))
        changed = True                             # every feature whose Parent chain leads to a removed gene
        while changed:
            changed = False
            for line in lines:
                f = line.rstrip("\n").split("\t")
                if len(f) > 8 and f[2] != "gene":
                    i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
                    if i and p and i.group(1) not in drop and set(p.group(1).split(",")) & drop:
                        drop.add(i.group(1))
                        changed = True
        with open(os.path.join(a.out_dir, "gff", os.path.basename(path)), "w") as fh:
            for line in lines:
                f = line.rstrip("\n").split("\t")
                if len(f) > 8:
                    i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
                    if (i and i.group(1) in drop) or (p and set(p.group(1).split(",")) & drop):
                        continue
                fh.write(line)

    # 2. matrix
    header, genomes, order, cells, extra = read_matrix(a.genes_matrix)
    gone = {g for _, g, _, _ in removed}
    pang_of = {}
    for p in order:
        for g in list(cells[p]):
            for x in cells[p][g]:
                pang_of[x] = p
            cells[p][g] = [x for x in cells[p][g] if x not in gone]
            if not cells[p][g]:
                del cells[p][g]
    emptied = [p for p in order if not cells[p]]
    order = [p for p in order if cells[p]]
    write_matrix(os.path.join(a.out_dir, "pangene_matrix_genes"), header, genomes, order, cells, extra)

    with open(os.path.join(a.out_dir, "removed_genes.tsv"), "w") as fh:
        fh.write("genome\tgene\tpanel_gene\tlocation\tpangene\tpangene_removed\treason\n")
        for genome, gid, b, loc in removed:
            p = pang_of.get(gid, "")
            fh.write(f"{genome}\t{gid}\t{b}\t{loc}\t{p or 'no pangene'}\t{'yes' if p in emptied else 'no'}\t"
                     f"{reason[b]}\n")
    if a.absences:
        with open(a.absences) as fi, open(os.path.join(a.out_dir, "absences_explained.tsv"), "w") as fo:
            for i, line in enumerate(fi):
                if i == 0 or line.split("\t", 1)[0] not in emptied:
                    fo.write(line)
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        fh.write(f"panel_genes_to_remove\t{len(reason)}\n"
                 f"panel_genes_found\t{len({b for _, _, b, _ in removed})}\ngenes_removed\t{len(removed)}\n"
                 f"genes_not_in_matrix\t{sum(1 for _, g, _, _ in removed if g not in pang_of)}\n"
                 f"pangenes_before\t{len(order) + len(emptied)}\npangenes_removed\t{len(emptied)}\n"
                 f"pangenes_after\t{len(order)}\n")
        for g, n in sorted(Counter(r[0] for r in removed).items()):
            fh.write(f"genes_removed_{g}\t{n}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read())


if __name__ == "__main__":
    main()
