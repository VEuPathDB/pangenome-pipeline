#!/usr/bin/env python3
"""
Find overlapping gene-level features in a final merged GFF3, for manual QC. Runs automatically
at the end of every genome's build (build_final_annotation.sh), so <genome>_overlaps.tsv is
ready to look at as soon as the build finishes - no need to wait for or run step 6's QC script
separately to get a first, per-pair look at what is overlapping and why.

agat_sp_complement_annotations.pl only rejects an --add gene when its CDS
overlaps the CDS of a same-strand gene already in --ref - gene-body/UTR overlap
alone is not treated as a conflict. Opposite-strand CDS overlaps (real antisense
pairs, but also occasional artifacts from complementing many sources together)
are always kept too. This script only reports gene pairs with an actual
CDS-to-CDS collision (UTR-only or gene-body-only overlap is not reported);
same-strand CDS overlaps are listed first since those are the ones that
shouldn't normally happen and are most worth a manual look.

Each pair's row also carries what is needed to judge it by hand without going back to the GFF:
each gene's description, CDS length, and whether it has InterPro/PFAM domain evidence (same
"has a domain" check the filter step and step 6's QC use), plus the overlap severity - the
overlap as a % of the SMALLER gene's own CDS (same rule/thresholds as this project's other
overlap-severity checks, so directly comparable to them).

Usage:
    python3 check_overlapping_genes.py --gff STRAIN_final_annotation.gff --out overlaps.tsv
"""

import argparse
import re
import sys
from collections import defaultdict
from urllib.parse import unquote

ID_RE = re.compile(r"ID=([^;]+)")
PARENT_RE = re.compile(r"Parent=([^;]+)")
DESCRIPTION_RE = re.compile(r"(?:^|;)description=([^;]+)")
PRODUCT_RE = re.compile(r"(?:^|;)product=([^;]+)")


def valid_orf(attrs):
    """True/False from Liftoff's own valid_ORF=/valid_ORFs= flag on the gene line - None if
    neither is present, which is normal (not a defect) for a gene that never went through
    Liftoff: the reference's own "main" genes, and every genome's own "original annotation"
    layer, are used as-is and never get this flag computed for them. Written as "NA" in the
    output, not blank/False, so that case isn't misread as "checked and found invalid"."""
    if "valid_ORF=True" in attrs:
        return True
    if "valid_ORF=False" in attrs:
        return False
    m = re.search(r"valid_ORFs=(\d+)", attrs)
    if m:
        return int(m.group(1)) >= 1
    return None


def parse_gff(gff_path, gene_feature_types):
    genes = {}  # gene_id -> dict
    mrna_to_gene = {}  # mrna_id -> gene_id
    cds_segments = defaultdict(list)  # gene_id -> list of (start, end)

    # First pass: genes and mRNA->gene links.
    with open(gff_path) as gff_file:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqname, source, feature, start, end, _, strand = fields[:7]
            match = ID_RE.search(fields[8])
            feat_id = match.group(1) if match else None

            if feature in gene_feature_types and feat_id:
                desc_m = DESCRIPTION_RE.search(fields[8]) or PRODUCT_RE.search(fields[8])
                genes[feat_id] = {
                    "id": feat_id, "source": source, "strand": strand, "seqname": seqname,
                    "start": int(start), "end": int(end),
                    "description": unquote(desc_m.group(1)) if desc_m else "",
                    "has_domain": False, "valid_ORF": valid_orf(fields[8]),
                }
            elif feature == "mRNA" and feat_id:
                parent = PARENT_RE.search(fields[8])
                if parent:
                    mrna_to_gene[feat_id] = parent.group(1)

    # Second pass: individual CDS segments, rolled up from mRNA to gene (kept
    # as separate segments, not one min/max span - a multi-exon gene's span
    # can cover an intron where another gene's CDS legitimately sits without
    # any real CDS-to-CDS collision, which a span-only check would miss).
    with open(gff_path) as gff_file:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            parent = PARENT_RE.search(fields[8])
            if not parent:
                continue
            gene_id = mrna_to_gene.get(parent.group(1))
            if gene_id not in genes:
                continue
            cds_segments[gene_id].append((int(fields[3]), int(fields[4])))
            if "InterPro:" in fields[8] or "PFAM:" in fields[8] or "Pfam:" in fields[8]:
                genes[gene_id]["has_domain"] = True

    return genes, cds_segments


def segments_overlap_bp(segments1, segments2):
    total = 0
    for s1, e1 in segments1:
        for s2, e2 in segments2:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi >= lo:
                total += hi - lo + 1
    return total


def cds_bp(gene_id, cds_segments):
    return sum(e - s + 1 for s, e in cds_segments.get(gene_id, []))


def overlap_severity(overlap_bp, cds_bp_a, cds_bp_b):
    """Same rule as this project's other overlap-severity checks (analyze_input_gffs.py,
    plot_final_annotation_qc.py): the overlap as a fraction of the SMALLER gene's own CDS
    length, not the bigger one and not the two combined."""
    smaller_bp = min(cds_bp_a, cds_bp_b)
    if smaller_bp == 0:
        return ""
    pct = overlap_bp / smaller_bp
    if pct >= 0.8:
        return "contained"
    elif pct >= 0.2:
        return "substantial"
    return "minor"


def find_overlaps(genes, cds_segments):
    by_seq = defaultdict(list)
    for gene in genes.values():
        by_seq[gene["seqname"]].append(gene)

    overlaps = []
    for seqname, seq_genes in by_seq.items():
        seq_genes.sort(key=lambda g: g["start"])
        active = []
        for gene in seq_genes:
            active = [a for a in active if a["end"] >= gene["start"]]
            for other in active:
                ov_bp = segments_overlap_bp(cds_segments.get(other["id"], []), cds_segments.get(gene["id"], []))
                overlaps.append({
                    "seqname": seqname,
                    "same_strand": other["strand"] == gene["strand"],
                    "overlap_bp": ov_bp,
                    "gene1": other, "gene2": gene,
                })
            active.append(gene)
    return overlaps


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True, help="Final merged GFF3 to check")
    parser.add_argument("--out", required=True, help="Output TSV, one row per overlapping gene pair")
    parser.add_argument("--gene-feature-types", default="gene,protein_coding_gene",
                         help="Comma-separated level-1 feature type name(s) that mark a gene (default: 'gene,protein_coding_gene')")
    args = parser.parse_args()

    gene_feature_types = {t.strip() for t in args.gene_feature_types.split(",") if t.strip()}

    genes, cds_segments = parse_gff(args.gff, gene_feature_types)
    overlaps = [ov for ov in find_overlaps(genes, cds_segments) if ov["overlap_bp"] > 0]
    overlaps.sort(key=lambda ov: not ov["same_strand"])

    columns = ["seqname", "strand_type", "severity", "overlap_bp",
               "gene1_id", "gene1_description", "gene1_source", "gene1_strand",
               "gene1_start", "gene1_end", "gene1_cds_bp", "gene1_has_domain", "gene1_valid_ORF",
               "gene2_id", "gene2_description", "gene2_source", "gene2_strand",
               "gene2_start", "gene2_end", "gene2_cds_bp", "gene2_has_domain", "gene2_valid_ORF"]

    with open(args.out, "w") as out_file:
        out_file.write("\t".join(columns) + "\n")
        for ov in overlaps:
            g1, g2 = ov["gene1"], ov["gene2"]
            bp1, bp2 = cds_bp(g1["id"], cds_segments), cds_bp(g2["id"], cds_segments)
            row = [
                ov["seqname"], "same_strand" if ov["same_strand"] else "opposite_strand",
                overlap_severity(ov["overlap_bp"], bp1, bp2), str(ov["overlap_bp"]),
                g1["id"], g1["description"], g1["source"], g1["strand"], str(g1["start"]), str(g1["end"]),
                str(bp1), str(g1["has_domain"]), "NA" if g1["valid_ORF"] is None else str(g1["valid_ORF"]),
                g2["id"], g2["description"], g2["source"], g2["strand"], str(g2["start"]), str(g2["end"]),
                str(bp2), str(g2["has_domain"]), "NA" if g2["valid_ORF"] is None else str(g2["valid_ORF"]),
            ]
            out_file.write("\t".join(row) + "\n")

    n_same = sum(1 for ov in overlaps if ov["same_strand"])
    print(f"[{args.gff}] {len(overlaps)} gene pairs with overlapping CDS "
          f"({n_same} same-strand - these shouldn't normally happen, worth a closer look) -> {args.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
