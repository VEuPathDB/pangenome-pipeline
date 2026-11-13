#!/usr/bin/env python3
"""
Build the per-genome gene index used by the explorer's locus view (the drawing of a pangene's genes and their
neighbours in every strain).

For every genome it writes <out>/loci_<genome>.json:
  {"contigs": {contig: [[gene_id, original_ID, source, start, end, strand, pangene, [cds_start, cds_end, ...], valid_orf,
                          excluded], ...]}}
with the genes of each contig sorted by start. pangene is the gene's pangene in --genes-matrix, or "" when the gene
is in the GFF but not in the matrix (e.g. removed by a filter). CDS segments are the union over the gene's
transcripts. valid_orf comes from step 05's orf_status: 1 (valid, repaired_miniprot, selenoprotein), 0 (invalid,
pseudogene) or -1 (ORF broken by the assembly, or no ORF information); without orf_status, from valid_ORFs.
excluded is the reason a gene was left out of the pangenome (its pangenome_excluded= tag, e.g. mito), or "".

Usage:
    python3 build_explorer_locus_index.py \
        --genes-matrix pangenes_final/matrix/pangene_matrix_genes.tr.tab \
        --gff-dir pangenes_final/gff --gff-suffix _final_annotation_tagged.gff -o pangenome_explorer/data

Needs only the standard library.
"""
import argparse
import csv
import json
import os
import re
import sys
import urllib.parse
from collections import defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
ORIG_RE = re.compile(r"(?:^|;)original_ID=([^;]+)")
VORF_RE = re.compile(r"(?:^|;)valid_ORFs=(\d+)")
STATUS_RE = re.compile(r"(?:^|;)orf_status=([^;]+)")
EXCL_RE = re.compile(r"(?:^|;)pangenome_excluded=([^;]+)")
ORF_CODE = {"valid": 1, "repaired_miniprot": 1, "selenoprotein": 1, "invalid": 0, "pseudogene": 0}
EMPTY = {"", "-", "NA"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("-o", "--out", required=True)
    args = ap.parse_args()

    with open(args.genes_matrix, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    cols = [(i, h.split("_for_get_pangenes")[0]) for i, h in enumerate(rows[0]) if i > 0 and h.strip()]
    pangene_of = defaultdict(dict)
    for r in rows[1:]:
        if not r or r[0].startswith("chr"):
            continue
        for i, g in cols:
            if i < len(r) and r[i] not in EMPTY:
                for x in r[i].split(","):
                    pangene_of[g][x] = r[0]

    os.makedirs(args.out, exist_ok=True)
    for _, g in cols:
        path = os.path.join(args.gff_dir, g + args.gff_suffix)
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")
        genes, parent, cds, ctg_of = {}, {}, defaultdict(set), {}
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
                if i and not p:
                    o, v, st = ORIG_RE.search(f[8]), VORF_RE.search(f[8]), STATUS_RE.search(f[8])
                    orf = (ORF_CODE.get(st.group(1).strip(), -1) if st else
                           (1 if int(v.group(1)) > 0 else 0) if v else -1)
                    ex = EXCL_RE.search(f[8])
                    genes[i.group(1)] = [i.group(1), o.group(1) if o else i.group(1), f[1], int(f[3]), int(f[4]), f[6],
                                         orf, urllib.parse.unquote(ex.group(1).strip()) if ex else ""]
                    ctg_of[i.group(1)] = f[0]
                elif i and p:
                    parent[i.group(1)] = p.group(1).split(",")[0]
                if f[2] == "CDS" and p:
                    cds[p.group(1).split(",")[0]].add((int(f[3]), int(f[4])))

        def root(x):
            seen = set()
            while x in parent and x not in seen:
                seen.add(x)
                x = parent[x]
            return x

        gene_cds = defaultdict(set)
        for tx, segs in cds.items():
            gene_cds[root(tx)] |= segs
        by_ctg = defaultdict(list)
        for gid, rec in genes.items():
            flat = [v for s, e in sorted(gene_cds.get(gid, ())) for v in (s, e)]
            by_ctg[ctg_of[gid]].append(rec[:6] + [pangene_of[g].get(gid, ""), flat, rec[6], rec[7]])
        for L in by_ctg.values():
            L.sort(key=lambda r: (r[3], r[4]))
        out = os.path.join(args.out, f"loci_{g}.json")
        with open(out, "w") as fo:
            json.dump({"contigs": by_ctg}, fo, separators=(",", ":"))
        n_out = sum(1 for L in by_ctg.values() for r in L if not r[6])
        n_ex = sum(1 for L in by_ctg.values() for r in L if r[9])
        print(f"{g}: {len(genes)} genes ({n_out} not in the matrix" + (f", {n_ex} of them left out of the pangenome" if n_ex else "")
              + f") -> {out} ({os.path.getsize(out) / 1e6:.1f} MB)", file=sys.stderr)


if __name__ == "__main__":
    main()
