#!/usr/bin/env python3
"""
Write an updated copy of each GFF3 with TE-derived gene calls reclassified
as standalone mobile_genetic_element features - the same shape RefSeq/GenBank
use for TE annotations (e.g. in Drosophila/human GFFs): column 3 (type) is
"mobile_genetic_element", with gbkey=mobile_element and a mobile_element_type
qualifier, and no mRNA/exon/CDS children (those aren't part of that
convention - a mobile_genetic_element is a flat feature, not a gene model).

Concretely, for each flagged gene: its gene/pseudogene line is replaced by
one mobile_genetic_element line at the same coordinates (fresh
ID=id-<contig>:<start>..<end>, gbkey=mobile_element,
mobile_element_type=<dominant_TE_class from te_gene_overlap.py>), and every
descendant line (mRNA/tRNA/rRNA, exon, CDS) is dropped. The original locus
tag, gene ID, gene_biotype, and TE_overlap_pct are kept as extra attributes
so the change is traceable back to the gene call it replaced. Every other
line in the file - every gene/feature not flagged - is copied through
unchanged.

Never touches --gff-indir: writes a complete copy of every genome's GFF
(edited where a gene was flagged, byte-identical where none were) to
--outdir, so --outdir ends up as a full, ready-to-use "updated annotation"
set alongside the untouched originals.

Which genes get reclassified:
  1. Automatically, every protein_coding gene with TE_overlap_pct >=
     --min-pct (default 80 - the "likely a TE" band, see the README).
  2. Optionally, genes named in --extra-genes: a 2-column TSV (no header)
     of genome_label<TAB>gene_id, for genes in the 20-80% "inspect" band
     that a manual look (e.g. at the extract_flagged_genes.py output)
     confirmed are TE-derived too. This is how a gene below --min-pct gets
     reclassified - nothing in the 20-80% band is touched automatically.

Usage:
    python3 -m pangenome_pipeline.step02_te_annotation_clean_up.reclassify_te_genes \
        --gene-table te_annotation_integration/gene_te_overlap.tsv \
        --gff-indir gff_originals \
        --outdir gff_updated \
        [--min-pct 80] [--extra-genes confirmed_20_80.tsv]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pangenome_pipeline.step02_te_annotation_clean_up.te_gene_overlap import (
    GENE_FTYPE_TO_BIOTYPE,
    parse_attributes,
    strip_prefix,
)

# Same set of gene-level feature types (GFF3 column 3) te_gene_overlap.py
# recognises when it builds gene_te_overlap.tsv in the first place - not just
# the NCBI-style "gene"/"pseudogene". VEuPathDB/FungiDB GFF3 (e.g. AF293 here)
# commonly uses "protein_coding_gene" instead of bare "gene"; a hardcoded
# ("gene", "pseudogene") tuple here would silently reclassify zero of its
# genes, leaving TE-derived gene models in place for exactly those genomes.
GENE_FTYPES = set(GENE_FTYPE_TO_BIOTYPE)

NEW_TYPE = "mobile_genetic_element"
FEATURE_TYPES_WITH_PARENT_GENE = {"mRNA", "tRNA", "rRNA"}
FEATURE_TYPES_WITH_PARENT_RNA = {"exon", "CDS"}


def load_extra_genes(path) -> set:
    """2-column TSV (no header): genome_label<TAB>gene_id."""
    pairs = set()
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genome_label, gene_id = line.split("\t")
            pairs.add((genome_label, gene_id))
    return pairs


def build_mobile_element_line(
    fields: list,
    gene_id: str,
    te_class: str,
    te_name: str,
    te_overlap_pct: float,
    original_type: str,
    original_attrs: dict,
) -> str:
    """fields: the original gene/pseudogene line's 9 columns (only 0
    (contig), 3 (start), 4 (end), 6 (strand) are reused - everything else is
    rebuilt fresh, matching how RefSeq/GenBank write mobile_genetic_element
    features (see e.g. Drosophila/human GFFs): no ID reused from the gene,
    a gbkey=mobile_element + mobile_element_type qualifier pair (formatted
    class:name, e.g. "transposon:1360" in the GenBank/Drosophila examples -
    here class is RepeatMasker's class/family, e.g. LTR/Gypsy, and name is
    the specific curated TE model, e.g. Gypsy-4_Pory), and no Parent/child
    structure."""
    contig, start, end, strand = fields[0], fields[3], fields[4], fields[6]
    locus_tag = original_attrs.get("locus_tag", gene_id)
    mobile_element_type = f"{te_class or 'unknown'}:{te_name or 'unknown'}"
    attrs = (
        f"ID=id-{contig}:{start}..{end}"
        f";gbkey=mobile_element"
        f";mobile_element_type={mobile_element_type}"
        f";original_locus_tag={locus_tag}"
        f";original_gene_id={gene_id}"
        f";original_type={original_type}"
        f";original_gene_biotype=protein_coding"
        f";TE_overlap_pct={te_overlap_pct:.2f}"
    )
    return "\t".join([contig, fields[1], NEW_TYPE, start, end, ".", strand, ".", attrs]) + "\n"


def reclassify_gff(gff_path: Path, gene_info: dict) -> tuple[list, int, set]:
    """gene_info: gene_id as it appears in gene_te_overlap.tsv (already
    stripped of the "gene-" prefix by te_gene_overlap.py) -> (TE_overlap_pct,
    dominant_TE_class), for every gene in this genome to reclassify.

    Two passes (order-independent w.r.t. parents appearing before
    children):
      1. Map every RNA feature's raw ID to its parent gene's raw ID (so an
         exon/CDS, whose Parent is an RNA ID, can be traced back to its
         gene), and find which genes' raw IDs are in gene_info.
      2. Walk the file again in original order: a flagged gene/pseudogene
         line is replaced by its mobile_genetic_element line; any
         mRNA/tRNA/rRNA/exon/CDS line belonging to a flagged gene is
         dropped; everything else is copied through unchanged.

    Returns (out_lines, n_changed, found_gene_ids) - found_gene_ids is the
    subset of gene_info's keys actually located in this GFF, so the caller
    can warn about the rest instead of silently leaving a TE-derived gene
    model in place."""
    lines = gff_path.read_text().splitlines(keepends=True)

    owning_gene_raw = {}
    target_raw = {}  # raw gene ID -> (gene_id, te_class, te_name, te_overlap_pct, original_type, attrs)
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
                gene_id = strip_prefix(fid, "gene-")
                if gene_id in gene_info:
                    pct, te_class, te_name = gene_info[gene_id]
                    target_raw[fid] = (gene_id, te_class, te_name, pct, ftype, attrs)
        elif ftype in FEATURE_TYPES_WITH_PARENT_GENE:
            parent = attrs.get("Parent")
            if fid and parent:
                owning_gene_raw[fid] = parent

    out_lines = []
    n_changed = 0
    for line in lines:
        if line.startswith("#"):
            out_lines.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            out_lines.append(line)
            continue
        ftype = fields[2]
        attrs = parse_attributes(fields[8])
        fid = attrs.get("ID")
        parent = attrs.get("Parent")

        if ftype in GENE_FTYPES:
            if fid in target_raw:
                gene_id, te_class, te_name, pct, original_type, original_attrs = target_raw[fid]
                out_lines.append(
                    build_mobile_element_line(fields, gene_id, te_class, te_name, pct, original_type, original_attrs)
                )
                n_changed += 1
            else:
                out_lines.append(line)
            continue

        if ftype in FEATURE_TYPES_WITH_PARENT_GENE:
            owner = parent
        elif ftype in FEATURE_TYPES_WITH_PARENT_RNA:
            owner = owning_gene_raw.get(parent)
        else:
            owner = None

        if owner in target_raw:
            continue  # drop: child of a reclassified gene
        out_lines.append(line)

    found_gene_ids = {gene_id for gene_id, *_ in target_raw.values()}
    return out_lines, n_changed, found_gene_ids


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-table", required=True, help="gene_te_overlap.tsv produced by te_gene_overlap.py")
    ap.add_argument(
        "--gff-indir",
        required=True,
        help="Directory with the source GFF3 files (genome_label.gff/.gff3) - never modified",
    )
    ap.add_argument(
        "--outdir",
        required=True,
        help="Directory to write the full updated GFF set (and a reclassified-genes summary TSV)",
    )
    ap.add_argument(
        "--min-pct",
        type=float,
        default=80.0,
        help="Reclassify every protein_coding gene with TE_overlap_pct >= this, automatically (default: 80)",
    )
    ap.add_argument(
        "--extra-genes",
        default=None,
        help="Optional 2-column TSV (genome_label<TAB>gene_id, no header) of "
        "manually-confirmed genes below --min-pct to reclassify as well",
    )
    args = ap.parse_args()

    gff_indir = Path(args.gff_indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.gene_table, sep="\t")
    auto = df[(df["gene_category"] == "protein_coding") & (df["TE_overlap_pct"] >= args.min_pct)]
    target = {
        (row.genome_label, row.gene_id): (row.TE_overlap_pct, row.dominant_TE_class, row.dominant_TE_name)
        for row in auto.itertuples()
    }

    if args.extra_genes:
        extra_pairs = load_extra_genes(args.extra_genes)
        indexed = df.set_index(["genome_label", "gene_id"])
        for key in extra_pairs:
            if key not in indexed.index:
                print(f"WARNING: --extra-genes entry {key} not found in {args.gene_table}, skipping")
                continue
            row = indexed.loc[key]
            target[key] = (row["TE_overlap_pct"], row["dominant_TE_class"], row["dominant_TE_name"])

    if not target:
        raise SystemExit(
            f"No genes to reclassify (TE_overlap_pct >= {args.min_pct}, plus none/no valid entries in --extra-genes)"
        )

    summary_rows = [
        {
            "genome_label": g,
            "gene_id": gid,
            "TE_overlap_pct": pct,
            "dominant_TE_class": te_class,
            "dominant_TE_name": te_name,
        }
        for (g, gid), (pct, te_class, te_name) in target.items()
    ]
    summary_df = pd.DataFrame(summary_rows).sort_values(["genome_label", "TE_overlap_pct"], ascending=[True, False])
    summary_path = outdir / "reclassified_genes.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"Wrote {summary_path} ({len(summary_df)} genes to reclassify)")

    by_genome_label = {}
    for (genome_label, gene_id), info in target.items():
        by_genome_label.setdefault(genome_label, {})[gene_id] = info

    gff_files = sorted(gff_indir.glob("*.gff")) + sorted(gff_indir.glob("*.gff3"))
    if not gff_files:
        raise SystemExit(f"ERROR: no .gff/.gff3 files found in {gff_indir}")

    total_changed = 0
    missing_pairs = []  # (genome_label, gene_id) targeted but never located/reclassified
    for gff_path in gff_files:
        genome_label = gff_path.stem
        gene_info = by_genome_label.get(genome_label, {})
        out_lines, n_changed, found_ids = reclassify_gff(gff_path, gene_info)
        missing_here = set(gene_info) - found_ids
        if missing_here:
            print(
                f"WARNING: {len(missing_here)} of {len(gene_info)} gene(s) not found in "
                f"{gff_path.name} for {genome_label!r} - NOT reclassified: {sorted(missing_here)}"
            )
            missing_pairs.extend((genome_label, gid) for gid in missing_here)
        dest = outdir / gff_path.name
        dest.write_text("".join(out_lines))
        total_changed += n_changed
        if n_changed:
            print(f"  {genome_label}: {n_changed} gene(s) reclassified -> {dest.name}")

    unmatched = set(by_genome_label) - {p.stem for p in gff_files}
    if unmatched:
        print(
            f"WARNING: {sum(len(by_genome_label[g]) for g in unmatched)} gene(s) in "
            f"genome_label(s) {sorted(unmatched)} had no matching GFF in {gff_indir} - not reclassified"
        )
        missing_pairs.extend((g, gid) for g in unmatched for gid in by_genome_label[g])

    if missing_pairs:
        print(
            f"\nWARNING: {len(missing_pairs)} of {len(target)} requested gene(s) were NOT "
            f"reclassified (still a gene model, not a mobile_genetic_element, in the output "
            f"GFF - see the WARNINGs above for why):"
        )
        for genome_label, gene_id in sorted(missing_pairs):
            print(f"  {genome_label}\t{gene_id}")

    print(
        f"\n{total_changed} gene(s) reclassified across {len(gff_files)} genome(s). "
        f"Full updated GFF set written to {outdir}/"
    )


if __name__ == "__main__":
    main()
