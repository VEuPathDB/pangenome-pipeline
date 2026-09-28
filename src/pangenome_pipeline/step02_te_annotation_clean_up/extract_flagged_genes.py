#!/usr/bin/env python3
"""
Extract the original GFF3 records for genes flagged by te_gene_overlap.py as
falling in a given TE-overlap band (by default, the 20-80% "inspect" band),
so they can be reviewed/curated directly instead of just as summary numbers.

For each matching genome, writes a subset GFF containing only the lines that
belong to the flagged genes - the gene/pseudogene feature itself plus every
descendant feature (mRNA/tRNA/rRNA, exon, CDS), in their original file order.
Nothing else from the source GFF is included.

Usage:
    python3 -m pangenome_pipeline.step02_te_annotation_clean_up.extract_flagged_genes \
        --gene-table te_annotation_integration/gene_te_overlap.tsv \
        --gff-indir gff_originals \
        --outdir te_annotation_integration/flagged_gff

Only extracts protein_coding and pseudogene genes - the two categories
te_gene_overlap.py scores (tRNA/rRNA genes are never in its table at all).
Defaults to 20% <= TE_overlap_pct < 80% (the "inspect" band from the
README); --min-pct/--max-pct override that - e.g. --min-pct 80 --max-pct
100.01 for the "likely a TE" band instead.

Expects --gene-table to be a gene_te_overlap.tsv produced by
te_gene_overlap.py (needs its genome_label and gene_id columns to locate
each gene's source GFF and record), and --gff-indir to contain those same
GFF files (genome_label.gff or .gff3).
"""

from __future__ import annotations

import argparse
from functools import partial
from pathlib import Path

import pandas as pd

from pangenome_pipeline.common.gff import find_gff as _find_gff
from pangenome_pipeline.step02_te_annotation_clean_up.te_gene_overlap import (
    GENE_FTYPE_TO_BIOTYPE,
    parse_attributes,
    strip_prefix,
)

# <genome_label>.gff/.gff3, falling back to a case-insensitive match for naming-convention mismatches.
find_gff = partial(_find_gff, exts=(".gff", ".gff3"), case_insensitive=True)

# The set of gene-level feature types (GFF3 column 3) recognised here must be
# the SAME set te_gene_overlap.py recognises when it builds gene_te_overlap.tsv
# in the first place (GENE_FTYPE_TO_BIOTYPE) - not just the NCBI-style "gene"/
# "pseudogene". VEuPathDB/FungiDB GFF3 (e.g. AF293 here) commonly uses
# "protein_coding_gene" instead of bare "gene"; a hardcoded ("gene",
# "pseudogene") tuple here silently recognised zero of its records, so every
# extracted GFF for such a genome came back empty even though the same genes
# were correctly scored in gene_te_overlap.tsv.
GENE_FTYPES = set(GENE_FTYPE_TO_BIOTYPE)

FEATURE_TYPES_WITH_PARENT_GENE = {"mRNA", "tRNA", "rRNA"}
FEATURE_TYPES_WITH_PARENT_RNA = {"exon", "CDS"}


ALWAYS_KEEP_HEADER_PREFIXES = ("##gff-version", "#!")


def extract_gene_records(gff_path: Path, target_gene_ids: set) -> tuple:
    """Two passes over the file's lines (order-independent w.r.t. parents
    appearing before children):
      1. Map every gene/pseudogene's raw ID to itself, and every RNA
         feature's raw ID to its parent gene's raw ID - together these let
         any exon/CDS (whose Parent is an RNA ID) be traced back to its gene.
         Also record which contigs those target genes actually sit on.
      2. Walk the file again in original order, keeping a feature line if
         the gene it belongs to is in target_gene_ids (matched after
         stripping the same "gene-" prefix te_gene_overlap.py strips when it
         builds gene_te_overlap.tsv, so the two stay consistent - the raw ID
         may or may not have had that prefix to begin with: NCBI GFFs do,
         VEuPathDB/FungiDB GFFs generally don't, and stripping a prefix that
         isn't there is a no-op, so this works for either). Header lines:
         the fixed top-of-file metadata (##gff-version, #!...) is always
         kept; a per-contig ##sequence-region/##species pair is only kept
         for contigs a target gene is actually on - these assemblies can
         have hundreds of scaffolds, and keeping every contig's header would
         swamp a handful of gene records in noise.

    Returns (out_lines, found_gene_ids) - found_gene_ids is the subset of
    target_gene_ids actually located in this GFF, so the caller can warn
    about the rest instead of just writing a silently-incomplete file."""
    lines = gff_path.read_text().splitlines(keepends=True)

    owning_gene_raw = {}
    target_raw_ids = set()
    target_contigs = set()
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        ftype = fields[2]
        attrs = parse_attributes(fields[8])
        fid = attrs.get("ID")
        if ftype in GENE_FTYPES:
            if fid:
                owning_gene_raw[fid] = fid
                if strip_prefix(fid, "gene-") in target_gene_ids:
                    target_raw_ids.add(fid)
                    target_contigs.add(fields[0])
        elif ftype in FEATURE_TYPES_WITH_PARENT_GENE:
            parent = attrs.get("Parent")
            if fid and parent:
                owning_gene_raw[fid] = parent

    out_lines = []
    current_seqregion_contig = None
    for line in lines:
        if line.startswith("#"):
            if line.startswith("##sequence-region"):
                parts = line.split()
                current_seqregion_contig = parts[1] if len(parts) > 1 else None
                if current_seqregion_contig in target_contigs:
                    out_lines.append(line)
            elif line.startswith("##species"):
                if current_seqregion_contig in target_contigs:
                    out_lines.append(line)
            elif line.startswith(ALWAYS_KEEP_HEADER_PREFIXES):
                out_lines.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        ftype = fields[2]
        attrs = parse_attributes(fields[8])
        fid = attrs.get("ID")
        parent = attrs.get("Parent")

        if ftype in GENE_FTYPES:
            owner = fid
        elif ftype in FEATURE_TYPES_WITH_PARENT_GENE:
            owner = parent
        elif ftype in FEATURE_TYPES_WITH_PARENT_RNA:
            owner = owning_gene_raw.get(parent)
        else:
            owner = None

        if owner in target_raw_ids:
            out_lines.append(line)
    found_gene_ids = {strip_prefix(fid, "gene-") for fid in target_raw_ids}
    return out_lines, found_gene_ids


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-table", required=True, help="gene_te_overlap.tsv produced by te_gene_overlap.py")
    ap.add_argument("--gff-indir", required=True, help="Directory with the source GFF3 files (genome_label.gff/.gff3)")
    ap.add_argument(
        "--outdir", required=True, help="Directory to write the flagged-genes summary TSV and per-genome GFF extracts"
    )
    ap.add_argument("--min-pct", type=float, default=20.0, help="Minimum TE_overlap_pct, inclusive (default: 20)")
    ap.add_argument("--max-pct", type=float, default=80.0, help="Maximum TE_overlap_pct, exclusive (default: 80)")
    ap.add_argument(
        "--min-div",
        type=float,
        default=None,
        help="Optional minimum TE_divergence_pct, inclusive - on top of --min-pct/--max-pct",
    )
    ap.add_argument("--max-div", type=float, default=None, help="Optional maximum TE_divergence_pct, exclusive")
    ap.add_argument(
        "--out-prefix",
        default="flagged_genes",
        help="Basename for the summary TSV (<prefix>.tsv) and per-genome GFF "
        "(<genome_label>.<prefix>.gff) - default: flagged_genes",
    )
    args = ap.parse_args()

    gff_indir = Path(args.gff_indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.gene_table, sep="\t")
    mask = (
        df["gene_category"].isin(("protein_coding", "pseudogene"))
        & (df["TE_overlap_pct"] >= args.min_pct)
        & (df["TE_overlap_pct"] < args.max_pct)
    )
    if args.min_div is not None:
        mask &= df["TE_divergence_pct"] >= args.min_div
    if args.max_div is not None:
        mask &= df["TE_divergence_pct"] < args.max_div
    sub = df[mask].sort_values(["genome_label", "TE_overlap_pct"], ascending=[True, False])

    if sub.empty:
        raise SystemExit(
            f"No protein_coding/pseudogene genes with {args.min_pct}% <= "
            f"TE_overlap_pct < {args.max_pct}% found in {args.gene_table}"
        )

    summary_path = outdir / f"{args.out_prefix}.tsv"
    sub.to_csv(summary_path, sep="\t", index=False)
    print(f"Wrote {summary_path} ({len(sub)} genes)")

    n_written = 0
    missing_pairs = []  # (genome_label, gene_id) never located in its GFF
    for genome_label, group in sub.groupby("genome_label"):
        gff_path = find_gff(gff_indir, genome_label)
        target_ids = set(group["gene_id"])
        if gff_path is None:
            print(
                f"WARNING: no GFF found for genome_label {genome_label!r} in {gff_indir} - skipping {len(group)} genes"
            )
            missing_pairs.extend((genome_label, gid) for gid in target_ids)
            continue
        out_lines, found_ids = extract_gene_records(gff_path, target_ids)
        missing_here = target_ids - found_ids
        if missing_here:
            print(
                f"WARNING: {len(missing_here)} of {len(target_ids)} gene(s) not found in "
                f"{gff_path.name} for {genome_label!r}: {sorted(missing_here)}"
            )
            missing_pairs.extend((genome_label, gid) for gid in missing_here)
        dest = outdir / f"{genome_label}.{args.out_prefix}.gff"
        dest.write_text("".join(out_lines))
        print(f"  {genome_label}: {len(found_ids)}/{len(target_ids)} genes -> {dest.name}")
        n_written += 1

    if missing_pairs:
        print(
            f"\nWARNING: {len(missing_pairs)} requested gene(s) across all genomes have "
            f"NO record in their output GFF (see the WARNINGs above for why):"
        )
        for genome_label, gene_id in sorted(missing_pairs):
            print(f"  {genome_label}\t{gene_id}")

    print(f"\nWrote GFF extracts for {n_written} genome(s) to {outdir}/")


if __name__ == "__main__":
    main()
