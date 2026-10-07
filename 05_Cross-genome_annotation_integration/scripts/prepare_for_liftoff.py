#!/usr/bin/env python3
"""
Prepare every genome's GFF for Liftoff. Two things, both done on a copy (the input GFFs are
never edited):

1. Choose what Liftoff transfers, by unifying the feature type of each top-level gene.
   Liftoff only transfers top-level features whose type is in its list - "gene" by default,
   plus whatever a `-f` types file adds (exact match, so protein_coding_gene, pseudogene or
   ncRNA_gene are NOT transferred unless listed). Annotations name their genes in different
   ways, so here every gene is classified (same rules as analyze_input_gffs.py: by its
   feature type, else gene_biotype=/biotype=/ebi_biotype=, else by what is nested under it)
   and retyped:
     - protein-coding genes and pseudogenes -> --gene-type (default protein_coding_gene) - what
       Liftoff will be told to transfer (see liftoff_feature_types.txt, written to --out-dir);
     - non-coding RNA genes (tRNA, rRNA, other ncRNA) are REMOVED from the prepared GFF, with
       every transcript and exon under them. The build step (AGAT) retypes a tRNA_gene as a plain
       "gene", which get_pangenes would then cluster, so they are left out instead.
   For the protein-coding genes only the type column changes; every attribute (including
   gene_biotype/ebi_biotype) is kept.
   Genes that can't be classified are left as they are and reported - if their type is
   literally "gene", Liftoff's default would still transfer them.

2. Remove the UTRs, so genes are described by their coding sequence only (nothing else gets
   transferred, and UTR overlap between neighbouring genes stops looking like gene overlap).
   The coding sequence is taken from the CDS lines themselves - not from any UTR line - so it
   works the same for annotations that list UTRs and ones that leave them implicit. For every
   transcript that has CDS lines (CDS extent = first to last CDS base):
     - UTR lines (five_prime_UTR, three_prime_UTR, ...) under it are deleted;
     - the transcript (mRNA) is shrunk to the CDS extent;
     - each exon is clipped to the CDS extent; an exon entirely outside it (UTR-only) is
       deleted; exons inside the CDS extent, and every CDS line, are left untouched;
     - each gene is shrunk to cover its transcripts.
   Genes/transcripts with no CDS (pseudogenes without a CDS) are not touched.

Usage:
    python3 prepare_for_liftoff.py \\
        --genomes genome_order.txt --gff-dir gff_curated --out-dir gff_for_liftoff

Writes <out-dir>/<genome>.gff3 (or .gff, as the input), liftoff_feature_types.txt (pass it to
Liftoff with -f) and prepare_for_liftoff_report.tsv.
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.gff import biotype_category, find_gff
from pangenome_common.io import read_id_list

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
BIOTYPE_RES = [
    re.compile(r"(?:^|;)gene_biotype=([^;]+)"),
    re.compile(r"(?:^|;)biotype=([^;]+)"),
    re.compile(r"(?:^|;)ebi_biotype=([^;]+)"),
]

TYPE_TO_CATEGORY = {
    "protein_coding_gene": "protein_coding",
    "pseudogene": "pseudogene",
    "trna_gene": "tRNA",
    "rrna_gene": "rRNA",
    "ncrna_gene": "ncRNA",
    "lncrna_gene": "ncRNA",
    "snrna_gene": "ncRNA",
    "snorna_gene": "ncRNA",
    "mirna_gene": "ncRNA",
    "scrna_gene": "ncRNA",
}
NCRNA_CHILD_TYPES = {"ncrna", "lnc_rna", "lncrna", "snrna", "snorna", "mirna", "misc_rna"}
NONCODING_GENE_TYPE = {"tRNA": "tRNA_gene", "rRNA": "rRNA_gene", "ncRNA": "ncRNA_gene"}


def is_utr(ftype):
    return "utr" in ftype.lower()


def classify(rec, children):
    t = rec["type"].lower()
    if t in TYPE_TO_CATEGORY:
        return TYPE_TO_CATEGORY[t]
    for rx in BIOTYPE_RES:
        m = rx.search(rec["f"][8])
        if m:
            return biotype_category(m.group(1))
    seen, stack, types = set(), list(children.get(rec["id"], [])), set()
    while stack:
        k = stack.pop()
        if id(k) in seen:
            continue
        seen.add(id(k))
        types.add(k["type"].lower())
        stack.extend(children.get(k["id"], []))
    if "cds" in types:
        return "protein_coding"
    if "trna" in types:
        return "tRNA"
    if "rrna" in types:
        return "rRNA"
    if types & NCRNA_CHILD_TYPES:
        return "ncRNA"
    if "pseudogenic_transcript" in types or "pseudogenic_exon" in types:
        return "pseudogene"
    return None


def prepare(gff_path, out_path, gene_type):
    with open(gff_path) as fh:
        lines = fh.readlines()

    recs, by_id = [], {}
    for line in lines:
        if line.startswith("#") or not line.strip():
            recs.append(None)
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            recs.append(None)
            continue
        m, pm = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        rec = {
            "f": f,
            "type": f[2],
            "start": int(f[3]),
            "end": int(f[4]),
            "id": m.group(1).strip() if m else None,
            "parents": [p.strip() for p in pm.group(1).split(",")] if pm else [],
            "drop": False,
        }
        recs.append(rec)
        if rec["id"]:
            by_id[rec["id"]] = rec

    children = defaultdict(list)
    for rec in recs:
        if rec:
            for p in rec["parents"]:
                children[p].append(rec)

    stats = defaultdict(int)
    retyped = defaultdict(int)  # (old type, new type) -> n
    unclassified = defaultdict(int)

    # ---- 1. unify gene types ----
    for rec in recs:
        if rec is None or rec["parents"] or "gene" not in rec["type"].lower():
            continue
        cat = classify(rec, children)
        old = rec["type"]
        if cat in ("protein_coding", "pseudogene"):
            new = gene_type
        elif cat in NONCODING_GENE_TYPE:
            # non-coding RNA genes are dropped, with everything nested under them
            stack = [rec]
            while stack:
                r = stack.pop()
                r["drop"] = True
                if r["id"]:
                    stack.extend(children.get(r["id"], []))
            retyped[(old, f"removed ({cat})")] += 1
            continue
        else:
            unclassified[old] += 1
            continue
        if new != old:
            rec["type"] = new
            retyped[(old, new)] += 1

    # ---- 2. remove UTRs ----
    cds_extent = {}
    for tid, kids in children.items():
        cds = [k for k in kids if k["type"] == "CDS"]
        if cds and tid in by_id:
            cds_extent[tid] = (min(k["start"] for k in cds), max(k["end"] for k in cds))

    for tid, (lo, hi) in cds_extent.items():
        by_id[tid]["start"], by_id[tid]["end"] = lo, hi
        stats["transcripts_trimmed"] += 1
        for kid in children[tid]:
            if is_utr(kid["type"]):
                kid["drop"] = True
                stats["utr_lines_deleted"] += 1
            elif kid["type"] == "exon":
                if not all(p in cds_extent for p in kid["parents"]):
                    continue
                elo = min(cds_extent[p][0] for p in kid["parents"])
                ehi = max(cds_extent[p][1] for p in kid["parents"])
                if kid["end"] < elo or kid["start"] > ehi:
                    if not kid["drop"]:
                        kid["drop"] = True
                        stats["exons_deleted"] += 1
                else:
                    ns, ne = max(kid["start"], elo), min(kid["end"], ehi)
                    if (ns, ne) != (kid["start"], kid["end"]) and not kid.get("clipped"):
                        kid["start"], kid["end"] = ns, ne
                        kid["clipped"] = True
                        stats["exons_clipped"] += 1

    for gid, kids in children.items():
        gene = by_id.get(gid)
        if gene is None or gid in cds_extent:
            continue
        if not any(k["id"] in cds_extent for k in kids):
            continue
        trs = [k for k in kids if k["type"] not in ("exon", "CDS") and not is_utr(k["type"])]
        new_lo, new_hi = min(k["start"] for k in trs), max(k["end"] for k in trs)
        if (new_lo, new_hi) != (gene["start"], gene["end"]):
            stats["genes_resized"] += 1
            stats["bp_trimmed_from_genes"] += (new_lo - gene["start"]) + (gene["end"] - new_hi)
            gene["start"], gene["end"] = new_lo, new_hi

    out = []
    for line, rec in zip(lines, recs):
        if rec is None:
            out.append(line)
        elif not rec["drop"]:
            f = rec["f"]
            f[2], f[3], f[4] = rec["type"], str(rec["start"]), str(rec["end"])
            out.append("\t".join(f) + "\n")
    with open(out_path, "w") as fh:
        fh.writelines(out)

    problems = 0
    for rec in recs:
        if rec and not rec["drop"]:
            for p in rec["parents"]:
                par = by_id.get(p)
                if (
                    par
                    and not par["drop"]
                    and rec["type"] in ("exon", "CDS")
                    and (rec["start"] < par["start"] or rec["end"] > par["end"])
                ):
                    problems += 1
    stats["problems"] = problems
    return stats, retyped, unclassified


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--genomes",
        required=True,
        help="Comma-separated genome names or a file with one per line (e.g. genome_order.txt)",
    )
    parser.add_argument("--gff-dir", required=True, help="Directory with each genome's <genome>.gff3 or <genome>.gff")
    parser.add_argument("--out-dir", required=True, help="Where to write the prepared GFFs (must not be --gff-dir)")
    parser.add_argument(
        "--gene-type",
        default="protein_coding_gene",
        help="Feature type given to every protein-coding gene and pseudogene, and the one Liftoff is "
        "told to transfer (default: protein_coding_gene - what the build/filter/QC scripts "
        "later in this step look for)",
    )
    args = parser.parse_args()

    if os.path.abspath(args.out_dir) == os.path.abspath(args.gff_dir):
        sys.exit("ERROR: --out-dir must be different from --gff-dir - the input GFFs are never edited in place")
    genomes = read_id_list(args.genomes)
    if not os.path.isdir(args.gff_dir):
        sys.exit(
            f"ERROR: --gff-dir {args.gff_dir!r} does not exist - did the previous step run, and are you in the project folder?"
        )
    missing = [g for g in genomes if find_gff(args.gff_dir, g) is None]
    if missing:
        sys.exit(
            f"ERROR: no GFF ({args.gff_dir}/<genome>.gff3 or .gff) for {len(missing)} of {len(genomes)} genome(s): "
            f"{', '.join(missing)}. Nothing was written - a genome silently missing here would be missing from everything after it."
        )

    os.makedirs(args.out_dir, exist_ok=True)

    with open(os.path.join(args.out_dir, "liftoff_feature_types.txt"), "w") as fh:
        fh.write(args.gene_type + "\n")

    cols = [
        "transcripts_trimmed",
        "utr_lines_deleted",
        "exons_clipped",
        "exons_deleted",
        "genes_resized",
        "bp_trimmed_from_genes",
    ]
    rows, all_retyped = [], defaultdict(int)
    for genome in genomes:
        gff_path = find_gff(args.gff_dir, genome)
        out_path = os.path.join(args.out_dir, f"{genome}.gff3")  # always .gff3, whatever the input's extension
        s, retyped, unclassified = prepare(gff_path, out_path, args.gene_type)
        rows.append((genome, [s[c] for c in cols]))
        print(
            f"[{genome}] {s['transcripts_trimmed']} transcripts trimmed, {s['utr_lines_deleted']} UTR lines "
            f"deleted, {s['exons_clipped']} exons clipped, {s['exons_deleted']} exons deleted, "
            f"{s['genes_resized']} genes resized ({s['bp_trimmed_from_genes']} bp) -> {out_path}",
            file=sys.stderr,
        )
        for (old, new), n in sorted(retyped.items()):
            print(f"    {n} x {old} -> {new}", file=sys.stderr)
            all_retyped[(genome, old, new)] += n
        for old, n in sorted(unclassified.items()):
            print(
                f"    WARNING: {n} x {old} could not be classified and was left as it is"
                + (' - typed "gene", so Liftoff\'s default would still transfer it' if old == "gene" else ""),
                file=sys.stderr,
            )
        if s["problems"]:
            print(
                f"    WARNING: {s['problems']} exon/CDS line(s) still stick out of their parent after "
                f"trimming - check {out_path}",
                file=sys.stderr,
            )

    with open(os.path.join(args.out_dir, "prepare_for_liftoff_report.tsv"), "w") as fh:
        fh.write("genome\t" + "\t".join(cols) + "\n")
        for genome, vals in rows:
            fh.write(genome + "\t" + "\t".join(str(v) for v in vals) + "\n")
        fh.write("\ngenome\told_gene_type\tnew_gene_type\tn_genes\n")
        for (genome, old, new), n in sorted(all_retyped.items()):
            fh.write(f"{genome}\t{old}\t{new}\t{n}\n")
    print(
        f"\nDone - prepared GFFs, liftoff_feature_types.txt and prepare_for_liftoff_report.tsv in {args.out_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
