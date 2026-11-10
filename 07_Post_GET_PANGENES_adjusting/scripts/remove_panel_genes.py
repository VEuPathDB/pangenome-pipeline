#!/usr/bin/env python3
"""
Remove genes found to be false during curation from the final GFFs and from the pangene matrix. Two lists:

--remove        whole panel genes (e.g. transposable elements): columns panel_gene, reason. A panel gene is a step 05
                basket: every gene tagged basket=<panel gene>, in every genome, is removed.
--remove-genes  single genes (optional): columns gene, reason. Only that gene is removed; the rest of its panel gene
                stays. The gene is given by its final ID (ID=) or by the genome's own gene ID (original_ID=, which does
                not change when step 05 is rerun). An ID found in no GFF is reported.

Each gene is removed with all its features (mRNA/transcript, exon, CDS...: every line whose Parent chain leads to the
gene). The other lines of the GFFs are written unchanged. In the matrix the removed genes leave their pangene; a pangene
left with no gene is removed. '#' lines of both lists are ignored; a list with only its header removes nothing.

Outputs (in --out-dir):
  gff/<genome><suffix>             the GFFs without the removed genes
  pangene_matrix_genes.tr.tab/.tab the matrix without them (pangenes in rows / genomes in rows)
  removed_genes.tsv                one row per gene removed: genome, gene, panel gene, location, pangene, the list that
                                   removed it (panel_gene or gene), reason
  genes_not_found.tsv              (with --remove-genes) the IDs of the gene list found in no GFF
  absences_explained.tsv           (with --absences) the same table without the rows of the pangenes removed
  summary.tsv                      counts

Usage:
    python3 remove_panel_genes.py --remove genes_to_remove.tsv --remove-genes strain_genes_to_remove.tsv \\
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
ORIG_RE = re.compile(r"(?:^|;)original_ID=([^;\n]+)")


def read_list(path, key):
    """{key column: reason} of a removal list ('#' lines ignored)."""
    return {r[key].strip(): r["reason"] for r in csv.DictReader(
        (l for l in open(path) if not l.startswith("#")), delimiter="\t") if r.get(key, "").strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--remove", required=True, help="TSV: panel_gene, reason (whole panel genes)")
    ap.add_argument("--remove-genes", default="", help="TSV: gene, reason (single genes: final ID or original_ID)")
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--absences", default="", help="absences_explained.tsv of fix_pangenome_safe.py (optional)")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    reason = read_list(a.remove, "panel_gene")
    gene_reason = read_list(a.remove_genes, "gene") if a.remove_genes else {}
    print(f"{len(reason)} panel genes and {len(gene_reason)} single genes to remove", file=sys.stderr)
    found = set()                                  # IDs of the gene list found in a GFF

    # 1. GFFs
    os.makedirs(os.path.join(a.out_dir, "gff"), exist_ok=True)
    removed = []                                   # (genome, gene, panel gene, location, list, reason)
    for path in sorted(glob.glob(os.path.join(a.gff_dir, "*" + a.gff_suffix))):
        genome = os.path.basename(path)[:-len(a.gff_suffix)]
        lines = open(path).readlines()
        drop = set()
        for line in lines:                         # the genes of the panel genes to remove
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene":
                b = BASKET_RE.search(f[8])
                gid = ID_RE.search(f[8]).group(1)
                o = ORIG_RE.search(f[8])
                hit = [x for x in (gid, o.group(1) if o else None) if x in gene_reason]
                found.update(hit)
                pg = b.group(1) if b else gid
                if pg in reason:
                    why = ("panel_gene", reason[pg])
                elif hit:
                    why = ("gene", gene_reason[hit[0]])
                else:
                    continue
                drop.add(gid)
                removed.append((genome, gid, pg, f"{f[0]}:{f[3]}-{f[4]}{f[6]}") + why)
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
    gone = {r[1] for r in removed}
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
        fh.write("genome\tgene\tpanel_gene\tlocation\tpangene\tpangene_removed\tlist\treason\n")
        for genome, gid, b, loc, lst, why in removed:
            p = pang_of.get(gid, "")
            fh.write(f"{genome}\t{gid}\t{b}\t{loc}\t{p or 'no pangene'}\t{'yes' if p in emptied else 'no'}\t"
                     f"{lst}\t{why}\n")
    missing = [g for g in gene_reason if g not in found]
    if a.remove_genes:
        with open(os.path.join(a.out_dir, "genes_not_found.tsv"), "w") as fh:
            fh.write("gene\treason\n")
            fh.writelines(f"{g}\t{gene_reason[g]}\n" for g in missing)
        if missing:
            print(f"WARNING: {len(missing)} IDs of {a.remove_genes} found in no GFF (genes_not_found.tsv)",
                  file=sys.stderr)
    if a.absences:
        with open(a.absences) as fi, open(os.path.join(a.out_dir, "absences_explained.tsv"), "w") as fo:
            for i, line in enumerate(fi):
                if i == 0 or line.split("\t", 1)[0] not in emptied:
                    fo.write(line)
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        fh.write(f"panel_genes_to_remove\t{len(reason)}\n"
                 f"panel_genes_found\t{len({r[2] for r in removed if r[4] == 'panel_gene'})}\n"
                 f"single_genes_to_remove\t{len(gene_reason)}\nsingle_genes_not_found\t{len(missing)}\n"
                 f"genes_removed\t{len(removed)}\n"
                 f"genes_removed_by_panel_gene_list\t{sum(r[4] == 'panel_gene' for r in removed)}\n"
                 f"genes_removed_by_gene_list\t{sum(r[4] == 'gene' for r in removed)}\n"
                 f"genes_not_in_matrix\t{sum(1 for r in removed if r[1] not in pang_of)}\n"
                 f"pangenes_before\t{len(order) + len(emptied)}\npangenes_removed\t{len(emptied)}\n"
                 f"pangenes_after\t{len(order)}\n")
        for g, n in sorted(Counter(r[0] for r in removed).items()):
            fh.write(f"genes_removed_{g}\t{n}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read())


if __name__ == "__main__":
    main()
