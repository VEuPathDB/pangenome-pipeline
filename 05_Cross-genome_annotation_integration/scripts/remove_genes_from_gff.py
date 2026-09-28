#!/usr/bin/env python3
"""
Remove specific genes (and every feature nested under them - mRNA, exon, CDS, UTR, ...) from
each genome's GFF, according to a manual kill list. Never edits the input GFFs: every genome
in --genomes gets a copy in --out-dir, with the listed genes stripped out of it (genomes with
nothing on the kill list are copied unchanged, so --out-dir is a complete drop-in replacement
for --gff-dir).

Kill list: TSV, no header, two columns - genome<TAB>gene_id. "#" lines and blank lines are
ignored. A gene ID only ever applies to the genome named on its own line, so the same ID under
two genomes needs two lines. Matches the feature's ID= attribute exactly; the NCBI-style
"gene-" and Ensembl-style "gene:" prefixed forms of the ID are tried too.

Usage:
    python3 remove_genes_from_gff.py \\
        --genomes genome_order.txt --gff-dir gff_with_interpro \\
        --kill-list kill_list_manual.tsv --out-dir gff_curated

Also writes <out-dir>/removed_genes_report.tsv: one row per removed gene, with its feature
type and how many GFF lines went with it. An ID on the kill list that isn't found in its
genome's GFF is reported as a WARNING (typo in the ID, or the wrong genome name) - nothing
else happens for it.
"""

import argparse
import os
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.gff import find_gff
from pangenome_common.io import read_id_list

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def read_kill_list(path):
    kills = defaultdict(set)
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or not parts[0].strip() or not parts[1].strip():
                sys.exit(f"ERROR: {path} line {n}: expected genome<TAB>gene_id, got {line!r}")
            kills[parts[0].strip()].add(parts[1].strip())
    return kills


def remove_genes(gff_path, out_path, gene_ids):
    with open(gff_path) as fh:
        lines = fh.readlines()

    id_to_type = {}
    children = defaultdict(list)
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        m = ID_RE.search(fields[8])
        if not m:
            continue
        feat_id = m.group(1).strip()
        id_to_type[feat_id] = fields[2]
        pm = PARENT_RE.search(fields[8])
        if pm:
            for parent in pm.group(1).split(","):
                children[parent.strip()].append(feat_id)

    removed_ids = set()
    matched = {}  # requested id -> actual ID found in the file
    for gid in gene_ids:
        for candidate in (gid, f"gene-{gid}", f"gene:{gid}"):
            if candidate in id_to_type:
                matched[gid] = candidate
                break
    for actual in matched.values():
        stack = [actual]
        while stack:
            fid = stack.pop()
            if fid in removed_ids:
                continue
            removed_ids.add(fid)
            stack.extend(children.get(fid, []))

    n_lines_per_gene = defaultdict(int)
    owner = {}
    for gid, actual in matched.items():
        stack = [actual]
        while stack:
            fid = stack.pop()
            if fid not in owner:
                owner[fid] = gid
                stack.extend(children.get(fid, []))

    kept = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            kept.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 9:
            m = ID_RE.search(fields[8])
            pm = PARENT_RE.search(fields[8])
            line_id = m.group(1).strip() if m else None
            parents = [p.strip() for p in pm.group(1).split(",")] if pm else []
            if line_id in removed_ids:
                n_lines_per_gene[owner[line_id]] += 1
                continue
            hit = next((p for p in parents if p in removed_ids), None)
            if hit is not None:
                n_lines_per_gene[owner[hit]] += 1
                continue
        kept.append(line)

    with open(out_path, "w") as fh:
        fh.writelines(kept)

    report = [(gid, id_to_type[matched[gid]], n_lines_per_gene[gid]) for gid in sorted(matched)]
    not_found = sorted(set(gene_ids) - set(matched))
    return report, not_found, len(lines) - len(kept)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genomes", required=True, help="Comma-separated genome names or a file with one per line (e.g. genome_order.txt)")
    parser.add_argument("--gff-dir", required=True, help="Directory with each genome's <genome>.gff3 or <genome>.gff")
    parser.add_argument("--kill-list", required=True, help="TSV, no header: genome<TAB>gene_id")
    parser.add_argument("--out-dir", required=True, help="Where to write the curated GFFs (must not be --gff-dir)")
    args = parser.parse_args()

    if os.path.abspath(args.out_dir) == os.path.abspath(args.gff_dir):
        sys.exit("ERROR: --out-dir must be different from --gff-dir - the input GFFs are never edited in place")

    genomes = read_id_list(args.genomes)
    if not os.path.isdir(args.gff_dir):
        sys.exit(f"ERROR: --gff-dir {args.gff_dir!r} does not exist - did the previous step run, and are you in the project folder?")
    missing = [g for g in genomes if find_gff(args.gff_dir, g) is None]
    if missing:
        sys.exit(f"ERROR: no GFF ({args.gff_dir}/<genome>.gff3 or .gff) for {len(missing)} of {len(genomes)} genome(s): "
                 f"{', '.join(missing)}. Nothing was written - a genome silently missing here would be missing from everything after it.")
    kills = read_kill_list(args.kill_list)
    unknown = sorted(set(kills) - set(genomes))
    if unknown:
        print(f"WARNING: kill list names genome(s) not in --genomes, ignored: {unknown}", file=sys.stderr)

    os.makedirs(args.out_dir, exist_ok=True)
    report_rows = []

    for genome in genomes:
        gff_path = find_gff(args.gff_dir, genome)
        out_path = os.path.join(args.out_dir, os.path.basename(gff_path))
        if genome not in kills:
            shutil.copyfile(gff_path, out_path)
            print(f"[{genome}] nothing on the kill list - copied unchanged", file=sys.stderr)
            continue
        report, not_found, n_removed = remove_genes(gff_path, out_path, kills[genome])
        for gid, ftype, n_lines in report:
            report_rows.append((genome, gid, ftype, n_lines))
        print(f"[{genome}] removed {len(report)} gene(s), {n_removed} GFF line(s) in total -> {out_path}", file=sys.stderr)
        for gid, ftype, n_lines in report:
            print(f"    {gid} ({ftype}): {n_lines} line(s)", file=sys.stderr)
        for gid in not_found:
            print(f"    WARNING: {gid} not found in {gff_path} - nothing removed for it (typo, or wrong genome?)", file=sys.stderr)

    with open(os.path.join(args.out_dir, "removed_genes_report.tsv"), "w") as fh:
        fh.write("genome\tgene_id\tfeature_type\tn_gff_lines_removed\n")
        for row in report_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"\nDone - curated GFFs + removed_genes_report.tsv in {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
