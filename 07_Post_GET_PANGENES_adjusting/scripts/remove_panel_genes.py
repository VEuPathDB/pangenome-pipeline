#!/usr/bin/env python3
"""
Remove genes found to be false during curation from the final GFFs and from the pangene matrix. Three lists:

--remove        whole panel genes (e.g. transposable elements): columns panel_gene, reason. A panel gene is a step 05
                basket: every gene tagged basket=<panel gene>, in every genome, is removed.
--remove-genes  single genes (optional): columns gene, reason. Only that gene is removed; the rest of its panel gene
                stays. The gene is given by its final ID (ID=) or by the genome's own gene ID (original_ID=, which does
                not change when step 05 is rerun). An ID found in no GFF is reported.
--exclude       panel genes left out of the pangenome only (optional; e.g. mitochondrial genes): columns panel_gene,
                reason. Every gene of the panel gene leaves the matrix and its counts, but stays in the GFFs (transfers
                and the strains' original genes alike), tagged pangenome_excluded=<reason>, so the explorer can still
                show it (in red). Which final GFF keeps them is decided later.
--actions       the table downloaded from the reference check page (optional): columns genome, item, target, panel_gene,
                action, reason. Its rows are added to the two lists above, and one more action:
                  delete                  -> A: the whole panel gene (panel_gene) is removed from every strain,
                                             original genes included, as in --remove
                  exclude_from_pangenome  -> B: every gene of the panel gene leaves the pangenome (matrix and counts);
                                             it stays in these GFFs, tagged pangenome_excluded=<reason>, as in --exclude
                  delete_in_genome        -> C: the gene (target) is not kept in its genome's reference annotation, but
                                             stays in the pangenome and in these GFFs, tagged reference_excluded=<reason>
                                             (listed in reference_excluded.tsv, for the official GFF)
                  keep_new_gene           -> nothing
                A wins over B and C; B wins over C.

Each gene is removed with all its features (mRNA/transcript, exon, CDS...: every line whose Parent chain leads to the
gene). The other lines of the GFFs are written unchanged. In the matrix the removed genes leave their pangene; a pangene
left with no gene is removed. '#' lines of both lists are ignored; a list with only its header removes nothing.
With --pangene-prefix the pangenes left are renamed <prefix>_pan_00001, <prefix>_pan_00002... in the order of the
matrix, and pangene_ids.tsv links each new ID to its old name (one strain's gene ID, given by get_pangenes).

Outputs (in --out-dir):
  gff/<genome><suffix>             the GFFs without the removed genes
  pangene_matrix_genes.tr.tab/.tab the matrix without them (pangenes in rows / genomes in rows), with the new pangene IDs
  pangene_ids.tsv                  (with --pangene-prefix) pangene_id, old_pangene, genomes, genes, panel_genes
  excluded_from_pangenome.tsv      (with --exclude or --actions) B: one row per gene left out of the pangenome: genome, gene, panel
                                   gene, location, pangene (old name), model_origin, reason
  reference_excluded.tsv           (with --actions) C: the genes not kept in their reference annotation
  removed_genes.tsv                one row per gene removed: genome, gene, panel gene, location, pangene (old name), the list
                                   that removed it (panel_gene or gene), reason
  genes_not_found.tsv              (with --remove-genes) the IDs of the gene list found in no GFF
  absences_explained.tsv           (with --absences) the same table without the rows of the pangenes removed, with the new
                                   pangene IDs (columns pangene, neighbours, match_pangene)
  summary.tsv                      counts

Usage:
    python3 remove_panel_genes.py --remove genes_to_remove.tsv --remove-genes strain_genes_to_remove.tsv \\
        --gff-dir basket_annotation/gff_prefer_panel_model_final \\
        --genes-matrix pangenome_fixed/pangene_matrix_genes.tr.tab \\
        --absences pangenome_fixed/absences_explained.tsv --pangene-prefix Tg --out-dir pangenome_final
"""
import argparse
import csv
import glob
import os
import re
import sys
import urllib.parse
from collections import Counter

from fix_pangenome_safe import read_matrix, write_matrix

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
BASKET_RE = re.compile(r"(?:^|;)basket=([^;\n]+)")
ORIG_RE = re.compile(r"(?:^|;)original_ID=([^;\n]+)")
MORIGIN_RE = re.compile(r"(?:^|;)model_origin=([^;\n]+)")


def read_list(path, key):
    """{key column: reason} of a removal list ('#' lines ignored)."""
    return {r[key].strip(): r["reason"] for r in csv.DictReader(
        (l for l in open(path) if not l.startswith("#")), delimiter="\t") if r.get(key, "").strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--remove", required=True, help="TSV: panel_gene, reason (whole panel genes)")
    ap.add_argument("--remove-genes", default="", help="TSV: gene, reason (single genes: final ID or original_ID)")
    ap.add_argument("--exclude", default="", help="TSV: panel_gene, reason (left out of the pangenome only)")
    ap.add_argument("--actions", default="", help="curation_actions.tsv of the reference check page (optional)")
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--absences", default="", help="absences_explained.tsv of fix_pangenome_safe.py (optional)")
    ap.add_argument("--pangene-prefix", default="", help="rename the pangenes <prefix>_pan_00001... (e.g. Tg, Af); "
                                                          "empty = keep get_pangenes' names")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    reason = read_list(a.remove, "panel_gene")
    gene_reason = read_list(a.remove_genes, "gene") if a.remove_genes else {}
    excl = {}                                      # B: panel genes left out of the pangenome
    ref_out = {}                                   # C: genes not kept in their reference annotation (kept here, tagged)
    if a.actions:
        for r in csv.DictReader((l for l in open(a.actions) if not l.startswith("#")), delimiter="\t"):
            act, why = r.get("action", "").strip(), (r.get("reason") or "").strip() or "curation"
            if act == "delete" and r.get("panel_gene"):
                reason.setdefault(r["panel_gene"].strip(), why)
            elif act == "delete_in_genome" and r.get("target"):
                ref_out.setdefault(r["target"].strip(), why)
            elif act == "exclude_from_pangenome" and r.get("panel_gene"):
                excl.setdefault(r["panel_gene"].strip(), why)
        ref_out = {k: v for k, v in ref_out.items() if k not in gene_reason}
    if a.exclude:
        for k, v in read_list(a.exclude, "panel_gene").items():
            excl.setdefault(k, v)
    excl = {k: v for k, v in excl.items() if k not in reason}
    print(f"{len(reason)} panel genes and {len(gene_reason)} single genes to remove; "
          f"{len(excl)} panel genes to leave out of the pangenome", file=sys.stderr)
    excluded = []                                  # (genome, gene, panel gene, location, reason, model_origin)
    basket_of = {}                                 # every gene -> its panel gene (for pangene_ids.tsv)
    refx = []                                      # (genome, gene, panel gene, location, reason)
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
                basket_of[gid] = pg
                if pg in reason:
                    why = ("panel_gene", reason[pg])
                elif hit:
                    why = ("gene", gene_reason[hit[0]])
                elif pg in excl:                   # out of the matrix only: stays in the GFF (tagged)
                    mo = MORIGIN_RE.search(f[8])
                    excluded.append((genome, gid, pg, f"{f[0]}:{f[3]}-{f[4]}{f[6]}", excl[pg], mo.group(1) if mo else ""))
                    continue
                else:
                    if gid in ref_out:
                        refx.append((genome, gid, pg, f"{f[0]}:{f[3]}-{f[4]}{f[6]}", ref_out[gid]))
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
        tag = {x[1]: "pangenome_excluded=" + urllib.parse.quote(x[4], safe=" ()/") for x in excluded if x[0] == genome}
        tag.update({x[1]: "reference_excluded=" + urllib.parse.quote(x[4], safe=" ()/") for x in refx if x[0] == genome})
        with open(os.path.join(a.out_dir, "gff", os.path.basename(path)), "w") as fh:
            for line in lines:
                f = line.rstrip("\n").split("\t")
                if len(f) > 8:
                    i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
                    if (i and i.group(1) in drop) or (p and set(p.group(1).split(",")) & drop):
                        continue
                    if f[2] == "gene" and i and i.group(1) in tag:   # left out of the pangenome: kept, tagged
                        f[8] = f[8].rstrip(";") + ";" + tag[i.group(1)]
                        line = "\t".join(f) + "\n"
                fh.write(line)

    # 2. matrix
    header, genomes, order, cells, extra = read_matrix(a.genes_matrix)
    gone = {r[1] for r in removed} | {x[1] for x in excluded}
    pang_of = {}
    for p in order:
        for g in list(cells[p]):
            for x in cells[p][g]:
                pang_of[x] = p
            cells[p][g] = [x for x in cells[p][g] if x not in gone]
            if not cells[p][g]:
                del cells[p][g]
    emptied = [p for p in order if not cells[p]]
    excl_genes = {x[1] for x in excluded}
    order = [p for p in order if cells[p]]
    new_id = {p: p for p in order}                 # old pangene name -> final ID
    if a.pangene_prefix:
        new_id = {p: f"{a.pangene_prefix}_pan_{i:05d}" for i, p in enumerate(order, 1)}
        with open(os.path.join(a.out_dir, "pangene_ids.tsv"), "w") as fh:
            fh.write("pangene_id\told_pangene\tgenomes\tgenes\tpanel_genes\n")
            for p in order:
                xs = [x for g in cells[p] for x in cells[p][g]]
                bs = sorted(Counter(basket_of.get(x, x) for x in xs).items(), key=lambda kv: (-kv[1], kv[0]))
                fh.write(f"{new_id[p]}\t{p}\t{len(cells[p])}\t{len(xs)}\t" + ",".join(b for b, _ in bs) + "\n")
        cells = {new_id[p]: cells[p] for p in order}
        order = [new_id[p] for p in order]
    write_matrix(os.path.join(a.out_dir, "pangene_matrix_genes"), header, genomes, order, cells, extra)

    with open(os.path.join(a.out_dir, "removed_genes.tsv"), "w") as fh:
        fh.write("genome\tgene\tpanel_gene\tlocation\tpangene\tpangene_removed\tlist\treason\n")
        for genome, gid, b, loc, lst, why in removed:
            p = pang_of.get(gid, "")
            fh.write(f"{genome}\t{gid}\t{b}\t{loc}\t{p or 'no pangene'}\t{'yes' if p in emptied else 'no'}\t"
                     f"{lst}\t{why}\n")
    if a.actions or a.exclude:
        with open(os.path.join(a.out_dir, "excluded_from_pangenome.tsv"), "w") as fh:
            fh.write("genome\tgene\tpanel_gene\tlocation\tpangene\tmodel_origin\treason\n")
            for genome, gid, b, loc, why, mo in excluded:
                fh.write(f"{genome}\t{gid}\t{b}\t{loc}\t{pang_of.get(gid, '') or 'no pangene'}\t{mo}\t{why}\n")
    if a.actions:
        with open(os.path.join(a.out_dir, "reference_excluded.tsv"), "w") as fh:
            fh.write("genome\tgene\tpanel_gene\tlocation\tpangene\treason\n")
            for genome, gid, b, loc, why in refx:
                fh.write(f"{genome}\t{gid}\t{b}\t{loc}\t{pang_of.get(gid, '') or 'no pangene'}\t{why}\n")
    missing = [g for g in gene_reason if g not in found]
    if a.remove_genes:
        with open(os.path.join(a.out_dir, "genes_not_found.tsv"), "w") as fh:
            fh.write("gene\treason\n")
            fh.writelines(f"{g}\t{gene_reason[g]}\n" for g in missing)
        if missing:
            print(f"WARNING: {len(missing)} IDs of {a.remove_genes} found in no GFF (genes_not_found.tsv)",
                  file=sys.stderr)
    if a.absences:
        ren = lambda v: "|".join(new_id.get(x, x) for x in v.split("|"))
        with open(a.absences) as fi, open(os.path.join(a.out_dir, "absences_explained.tsv"), "w") as fo:
            for i, line in enumerate(fi):
                f = line.rstrip("\n").split("\t")
                if i == 0:
                    cols = [j for j, h in enumerate(f) if h in ("pangene", "neighbours", "match_pangene")]
                elif f[0] in emptied:
                    continue
                else:
                    for j in cols:
                        if j < len(f) and f[j]:
                            f[j] = ren(f[j])
                fo.write("\t".join(f) + "\n")
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        fh.write(f"panel_genes_to_remove\t{len(reason)}\n"
                 f"panel_genes_found\t{len({r[2] for r in removed if r[4] == 'panel_gene'})}\n"
                 f"single_genes_to_remove\t{len(gene_reason)}\nsingle_genes_not_found\t{len(missing)}\n"
                 f"genes_removed\t{len(removed)}\n"
                 f"genes_removed_by_panel_gene_list\t{sum(r[4] == 'panel_gene' for r in removed)}\n"
                 f"genes_removed_by_gene_list\t{sum(r[4] == 'gene' for r in removed)}\n"
                 f"panel_genes_left_out_of_pangenome\t{len({x[2] for x in excluded})}\n"
                 f"genes_left_out_of_pangenome\t{len(excluded)}\n"
                 f"genes_left_out_original\t{sum(1 for x in excluded if x[5] == 'original')}\n"
                 f"genes_not_kept_in_reference\t{len(refx)}\n"
                 f"genes_not_in_matrix\t{sum(1 for r in removed if r[1] not in pang_of)}\n"
                 f"pangenes_before\t{len(order) + len(emptied)}\npangenes_removed\t{len(emptied)}\n"
                 f"pangenes_left_out\t{sum(1 for p in emptied if any(pang_of.get(g) == p for g in excl_genes))}\n"
                 f"pangenes_after\t{len(order)}\n")
        for g, n in sorted(Counter(r[0] for r in removed).items()):
            fh.write(f"genes_removed_{g}\t{n}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read())


if __name__ == "__main__":
    main()
