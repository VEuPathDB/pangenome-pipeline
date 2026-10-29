#!/usr/bin/env python3
"""
The fate of every native gene in the final GFFs: the final model(s) at its place, or the record that explains why it
has none. No gene may be lost without a record.

A native gene (ID, no Parent, with a CDS) is "kept" when a final gene of its own basket (a model of one of its basket's
winners) overlaps its CDS span. A model of another panel gene at its place does not keep it: it is listed in
overlapped_by. Otherwise its reason is looked up, in this order:
  left out                left_out.tsv (build): its basket has a winner, no model of it covers the gene, and the gene
                          is a piece of the basket's model or failed the check
  not eligible            not_eligible.tsv (baskets_v2.py): could not win, and no other member of its basket won
  panel gene removed      removed_panel_genes.tsv (final): its basket's panel gene was removed (curation, overlap)
  one place, one gene     stacked_dropped.tsv (build): its basket's model was dropped for another panel gene's model
                          and the gene was not kept at another place
  no record               none of the above: an error (exit 1 with --strict)

Writes <out>: genome, gene, location, fate, detail (the model(s) that keep it, or the reason), overlapped_by (final
genes of other panel genes over its CDS span, either strand). Prints the counts per genome.

Usage:
    python3 gene_fate.py --native-dir basket_annotation/input/native --genomes genome_order.txt \
        --result-dir basket_annotation --build-dir basket_annotation/gff_prefer_panel_model \
        --final-dir basket_annotation/gff_prefer_panel_model_final --out basket_annotation/gene_fate.tsv --strict
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_basket_annotation import gene_blocks, span   # noqa: E402

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
BASKET_RE = re.compile(r"(?:^|;)basket=([^;\n]+)")


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--native-dir", required=True, help="the native GFFs given to the build (<genome>.gff3 / .gff)")
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--result-dir", required=True, help="baskets_v2.tsv, not_eligible.tsv")
    ap.add_argument("--build-dir", required=True, help="the build's output: left_out.tsv, stacked_dropped.tsv")
    ap.add_argument("--final-dir", required=True, help="final GFFs and removed_panel_genes.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--strict", action="store_true", help="exit 1 when a gene has no record")
    a = ap.parse_args()
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]

    winners_of = {}
    for r in rows(os.path.join(a.result_dir, "baskets_v2.tsv")):
        for m in r["members"].split(","):
            winners_of[m] = [w for w in r["winners"].split(";") if w]
    ne = {(r["genome"], r["gene"]): r["reason"] for r in rows(os.path.join(a.result_dir, "not_eligible.tsv"))}
    left = {(r["genome"], r["gene"]): r["reason"] for r in rows(os.path.join(a.build_dir, "left_out.tsv"))}
    removed = {r["panel_gene"]: r["reason"] for r in rows(os.path.join(a.final_dir, "removed_panel_genes.tsv"))}
    stacked = {(r["genome"], r["panel_gene_dropped"]): r["kept_panel_gene"]
               for r in rows(os.path.join(a.build_dir, "stacked_dropped.tsv"))}

    out, counts = [], defaultdict(Counter)
    for g in genomes:
        models = defaultdict(list)
        for line in open(os.path.join(a.final_dir, f"{g}_basket_annotation.gff")):
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                b = BASKET_RE.search(f[8])
                models[f[0]].append((int(f[3]), int(f[4]), ID_RE.search(f[8]).group(1), b.group(1) if b else "-"))
        path = next(iter(glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff"))))
        ids = set()
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and "Parent=" not in f[8]:
                m = ID_RE.search(f[8])
                if m:
                    ids.add(m.group(1))
        for gid, block in sorted(gene_blocks(path, ids).items()):
            c, st, s, e = span(block)
            ws = winners_of.get(gid, [])
            over = [m for m in models[c] if m[0] <= e and m[1] >= s]
            hit = [m for m in over if m[3] in ws]
            other = ";".join(f"{m[2]} ({m[3]})" for m in over if m[3] not in ws) or "-"
            if hit:
                fate, detail = "kept", ";".join(f"{m[2]} ({m[3]})" for m in hit)
            elif (g, gid) in left:
                fate, detail = "left out", left[(g, gid)]
            elif (g, gid) in ne:
                fate, detail = "not eligible", ne[(g, gid)]
            elif any(w in removed for w in ws):
                fate, detail = "panel gene removed", "; ".join(f"{w}: {removed[w]}" for w in ws if w in removed)
            elif any((g, w) in stacked for w in ws):
                fate, detail = "one place, one gene", "; ".join(f"{w} dropped, kept {stacked[(g, w)]}"
                                                                for w in ws if (g, w) in stacked)
            else:
                fate, detail = "NO RECORD", "basket " + (";".join(ws) or "without a winner")
            counts[g][fate] += 1
            out.append([g, gid, f"{c}:{s}-{e}({st})", fate, detail, other])

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "gene", "location", "fate", "detail", "overlapped_by"])
        w.writerows(out)
    fates = ["kept", "not eligible", "left out", "panel gene removed", "one place, one gene", "NO RECORD"]
    print("genome\t" + "\t".join(fates))
    for g in genomes:
        print(g + "\t" + "\t".join(str(counts[g][x]) for x in fates))
    lost = sum(c["NO RECORD"] for c in counts.values())
    print(f"genes without a record: {lost} ({a.out})")
    if lost and a.strict:
        sys.exit(f"ERROR: {lost} native genes have no final model and no record (fate NO RECORD in {a.out})")


if __name__ == "__main__":
    main()
