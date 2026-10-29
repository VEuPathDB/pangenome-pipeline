#!/usr/bin/env python3
"""
Pilot: add the VEuPathDB OrthoMCL orthogroups of the reference genes to the explorer data (nothing else is changed).

A pangene gets the orthogroup of every reference gene it holds: its panel genes (step 05 basket=) and the original IDs
of its reference-strain genes. Pangenes with no reference gene (panel gene from another strain, new genes) get none.

Input: the "Genes by Ortholog Pattern" download of the reference from ToxoDB/FungiDB (tab-separated, columns
"Gene ID" and "Ortholog Group").

Output (--out-dir): a copy of the explorer folder where every pangene has "orthogroups": [[group, pangenes in it]],
and orthogroups_per_pangene.tsv / summary.tsv.

Usage:
    python3 add_orthogroups.py --explorer-dir pangenome_explorer/preliminary \\
        --og-table GenesByOrthologPattern_Summary.txt --reference ME49 \\
        --out-dir pangenome_explorer/preliminary_orthogroups
"""
import argparse
import csv
import glob
import json
import os
import shutil
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--explorer-dir", required=True, help="explorer folder with data/ (build_explorer_job.sh output)")
    ap.add_argument("--og-table", required=True, help="VEuPathDB Genes by Ortholog Pattern download of the reference")
    ap.add_argument("--reference", required=True, help="reference strain name in the explorer (ME49, AF293)")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    og = {r["Gene ID"]: r["Ortholog Group"] for r in csv.DictReader(open(a.og_table), delimiter="\t")
          if r.get("Ortholog Group")}
    src = os.path.join(a.explorer_dir, "data")
    man = json.load(open(os.path.join(src, "manifest.json")))
    chunks = sorted(glob.glob(os.path.join(src, "chunk_*.json")))
    data = [json.load(open(f)) for f in chunks]

    # field positions of a gene record (same order build_explorer_data.py writes)
    fields = man.get("fields") or man.get("gene_fields")
    oi = fields.index("original_ID") if fields and "original_ID" in fields else None
    bi = fields.index("basket") if fields and "basket" in fields else None

    by_og, per = defaultdict(set), {}
    for chunk in data:
        for p in chunk:
            ids = {x[0] for x in p.get("panel_genes", [])}
            for rec in p["genes"].get(a.reference, []):
                for i in (oi, bi):
                    if i is not None and rec[i]:
                        ids.add(rec[i])
            groups = sorted({og[i] for i in ids if i in og})
            per[p["id"]] = groups
            for g in groups:
                by_og[g].add(p["id"])

    if os.path.exists(a.out_dir):
        shutil.rmtree(a.out_dir)
    shutil.copytree(a.explorer_dir, a.out_dir, ignore=shutil.ignore_patterns("explorer_standalone*.html", "chunk_*.json"))
    for f, chunk in zip(chunks, data):
        for p in chunk:
            p["orthogroups"] = [[g, len(by_og[g])] for g in per[p["id"]]]
        with open(os.path.join(a.out_dir, "data", os.path.basename(f)), "w") as fh:
            json.dump(chunk, fh, separators=(",", ":"))

    with open(os.path.join(a.out_dir, "orthogroups_per_pangene.tsv"), "w") as fh:
        fh.write("pangene\torthogroups\tpangenes_in_orthogroup\n")
        for pid, groups in per.items():
            fh.write(f"{pid}\t{';'.join(groups)}\t{';'.join(str(len(by_og[g])) for g in groups)}\n")
    n = len(per)
    with_og = sum(1 for g in per.values() if g)
    summ = [("pangenes", n), ("pangenes_with_orthogroup", with_og), ("pangenes_without_orthogroup", n - with_og),
            ("pangenes_with_2+_orthogroups", sum(1 for g in per.values() if len(g) > 1)),
            ("orthogroups", len(by_og)),
            ("orthogroups_in_2+_pangenes", sum(1 for v in by_og.values() if len(v) > 1)),
            ("pangenes_sharing_an_orthogroup", len({p for v in by_og.values() if len(v) > 1 for p in v}))]
    with open(os.path.join(a.out_dir, "orthogroups_summary.tsv"), "w") as fh:
        for k, v in summ:
            fh.write(f"{k}\t{v}\n")
            print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
