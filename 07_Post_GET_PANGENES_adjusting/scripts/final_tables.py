#!/usr/bin/env python3
"""
Two tables of the final pangenome (after remove_panel_genes.py):

  gene_status.tsv   every gene of every strain in step 05's GFFs, and what happened to it after step 05: its pangene in
                    get_pangenes and in the final pangenome, and its status:
                      in_pangenome         in the same pangene as in get_pangenes
                      moved                moved to another pangene by fix_pangenome_safe.py (same gene, split by
                                           get_pangenes)
                      not_in_pangenome     in the final GFF but in no pangene (get_pangenes leaves out CDS < 100 bp)
                      removed_curation     removed with its whole panel gene by remove_panel_genes.py
                                           (genes_to_remove.tsv; reason given)
                      removed_curation_gene  removed alone by remove_panel_genes.py (strain_genes_to_remove.tsv;
                                           reason given)
                      excluded_from_pangenome  left out of the pangenome with its whole panel gene
                                           (genes_to_exclude_from_pangenome.tsv or curation_actions.tsv,
                                           exclude_from_pangenome); still in the final GFFs; note: its model origin and
                                           the reason
                      not_in_reference     in the pangenome and in these GFFs, but not kept in its reference annotation
                                           (curation_actions.tsv, delete_in_genome; reason given)
                    plus one row per panel gene removed or left out by step 05 itself (genome "-": removed_step05,
                    not_eligible_step05), so every gene of the panel is accounted for.
  singletons.tsv    every pangene with a gene in a single strain, and why it is alone (first match):
                      tandem_copy          a tandem-array copy added by step 05
                      split_by_get_pangenes  its panel gene also has genes in other pangenes (get_pangenes separated it)
                      inside_another_gene  its span lies >= 90% inside another gene of the same strain (either strand)
                      partial              its CDS is < 50% of its panel gene's protein length
                      unique_to_strain     its panel gene has no gene in any other strain
                      other                none of the above
  summary.tsv       counts

Usage:
    python3 final_tables.py --step05-gff-dir basket_annotation/gff_prefer_panel_model_final \\
        --step05-dir basket_annotation --pangenes-matrix <get_pangenes>/pangene_matrix_genes.tr.tab \\
        --moves pangenome_fixed/moves.tsv --final-dir pangenome_final
"""
import argparse
import csv
import os
from collections import Counter, defaultdict

from fix_pangenome_safe import read_matrix
from summarize_pangenes import read_genes


def gene_to_pangene(path):
    _, _, order, cells, _ = read_matrix(path)
    return {x: p for p in order for xs in cells[p].values() for x in xs}, cells


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step05-gff-dir", required=True, help="step 05 final GFFs (before step 07)")
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--step05-dir", required=True, help="step 05 basket_annotation folder (winners.tsv, not_eligible.tsv, "
                                                        "gff_prefer_panel_model_final/removed_panel_genes.tsv)")
    ap.add_argument("--pangenes-matrix", required=True, help="get_pangenes' pangene_matrix_genes.tr.tab")
    ap.add_argument("--moves", required=True, help="fix_pangenome_safe.py's moves.tsv")
    ap.add_argument("--final-dir", required=True, help="remove_panel_genes.py's output (matrix, removed_genes.tsv)")
    a = ap.parse_args()

    genes = read_genes(a.step05_gff_dir, a.gff_suffix)
    tandem = set()
    for path in [os.path.join(a.step05_gff_dir, f) for f in os.listdir(a.step05_gff_dir) if f.endswith(a.gff_suffix)]:
        for line in open(path):
            if "\tgene\t" in line and "tandem_array_copy=True" in line:
                tandem.add(line.split("ID=", 1)[1].split(";", 1)[0].strip())
    win = {r["gene"]: r for r in csv.DictReader(open(os.path.join(a.step05_dir, "winners.tsv")), delimiter="\t")}
    first, _ = gene_to_pangene(a.pangenes_matrix)
    final, cells = gene_to_pangene(os.path.join(a.final_dir, "pangene_matrix_genes.tr.tab"))
    moved = {r["gene"]: r["from_pangene"] for r in csv.DictReader(open(a.moves), delimiter="\t")}
    removed = {r["gene"]: (r.get("list", "panel_gene"), r["reason"]) for r in
               csv.DictReader(open(os.path.join(a.final_dir, "removed_genes.tsv")), delimiter="\t")}
    ex_path = os.path.join(a.final_dir, "excluded_from_pangenome.tsv")
    excluded = {r["gene"]: (r.get("model_origin") or "-") + "; " + r["reason"]
                for r in csv.DictReader(open(ex_path), delimiter="\t")} if os.path.exists(ex_path) else {}
    rx_path = os.path.join(a.final_dir, "reference_excluded.tsv")
    ref_out = {r["gene"]: r["reason"] for r in csv.DictReader(open(rx_path), delimiter="\t")} if os.path.exists(rx_path) else {}

    # 1. every gene
    rows = []
    for g, x in sorted(genes.items(), key=lambda kv: (kv[1]["genome"], kv[0])):
        if g in excluded:
            status, note = "excluded_from_pangenome", excluded[g]
        elif g in removed:
            status = "removed_curation" if removed[g][0] == "panel_gene" else "removed_curation_gene"
            note = removed[g][1]
        elif g in ref_out:
            status, note = "not_in_reference", ref_out[g]
        elif g in final:
            status, note = ("moved", f"from {moved[g]}") if g in moved else ("in_pangenome", "")
        else:
            status, note = "not_in_pangenome", f"CDS {x['cds_bp']} bp"
        rows.append([x["genome"], g, x["basket"], x["location"], x["model_origin"], x["orf_status"], x["cds_bp"],
                     win.get(x["basket"], {}).get("psauron", ""), first.get(g, ""), final.get(g, ""), status, note])
    rp = os.path.join(a.step05_dir, "gff_prefer_panel_model_final", "removed_panel_genes.tsv")
    if os.path.exists(rp):
        for r in csv.DictReader(open(rp), delimiter="\t"):
            rows.append(["-", "-", r["panel_gene"], "", "", "", "", "", "", "", "removed_step05", r["reason"]])
    ne = os.path.join(a.step05_dir, "not_eligible.tsv")
    if os.path.exists(ne):
        for r in csv.DictReader(open(ne), delimiter="\t"):
            rows.append(["-", "-", r["gene"], "", "", "", "", r.get("psauron", ""), "", "", "not_eligible_step05",
                         r["reason"]])
    with open(os.path.join(a.final_dir, "gene_status.tsv"), "w") as fh:
        fh.write("genome\tgene\tpanel_gene\tlocation\tmodel_origin\torf_status\tcds_bp\tpsauron\t"
                 "pangene_get_pangenes\tpangene_final\tstatus\tnote\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    # 2. singletons
    by_basket = defaultdict(set)                   # panel gene -> pangenes / strains with its genes in the final pangenome
    strains_of = defaultdict(set)
    for g, p in final.items():
        by_basket[genes[g]["basket"]].add(p)
        strains_of[genes[g]["basket"]].add(genes[g]["genome"])
    spans = defaultdict(list)                      # genome, contig -> (start, end, gene)
    for g, x in genes.items():
        if g in removed:
            continue
        c, r = x["location"].rsplit(":", 1)
        s, e = r[:-1].split("-")
        spans[(x["genome"], c)].append((int(s), int(e), g))
    sing = []
    for p, cg in cells.items():
        if len(cg) != 1:
            continue
        for g in next(iter(cg.values())):
            x = genes[g]
            c, r = x["location"].rsplit(":", 1)
            s, e = map(int, r[:-1].split("-"))
            host = [h for hs, he, h in spans[(x["genome"], c)]
                    if h != g and min(e, he) - max(s, hs) + 1 >= 0.9 * (e - s + 1)]
            aa = win.get(x["basket"], {}).get("protein_aa", "")
            others = sorted(by_basket[x["basket"]] - {p})
            if g in tandem:
                cls, det = "tandem_copy", ""
            elif others:
                cls, det = "split_by_get_pangenes", ",".join(others)
            elif host:
                cls, det = "inside_another_gene", ",".join(f"{h} ({final.get(h, 'no pangene')})" for h in host)
            elif aa and x["cds_bp"] / 3 < 0.5 * int(aa):
                cls, det = "partial", f"CDS {x['cds_bp']} bp, panel gene {aa} aa"
            elif strains_of[x["basket"]] == {x["genome"]}:
                cls, det = "unique_to_strain", ""
            else:
                cls, det = "other", ""
            sing.append([p, x["genome"], g, x["basket"], cls, det, x["cds_bp"], x["orf_status"],
                         win.get(x["basket"], {}).get("psauron", ""), x["description"]])
    with open(os.path.join(a.final_dir, "singletons.tsv"), "w") as fh:
        fh.write("pangene\tgenome\tgene\tpanel_gene\tclass\tdetail\tcds_bp\torf_status\tpsauron\tdescription\n")
        for r in sorted(sing, key=lambda r: (r[4], r[1], r[0])):
            fh.write("\t".join(map(str, r)) + "\n")

    st, sc = Counter(r[10] for r in rows), Counter(r[4] for r in sing)
    with open(os.path.join(a.final_dir, "summary_tables.tsv"), "w") as fh:
        for k, n in st.most_common():
            fh.write(f"{'panel_genes' if k.endswith('_step05') else 'genes'}_{k}\t{n}\n")
        fh.write(f"singleton_pangenes\t{len({r[0] for r in sing})}\n")
        for k, n in sc.most_common():
            fh.write(f"singletons_{k}\t{n}\n")
    print(open(os.path.join(a.final_dir, "summary_tables.tsv")).read())


if __name__ == "__main__":
    main()
