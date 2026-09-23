#!/usr/bin/env python3
"""
Profile every genome's own input GFF before Liftoff runs - what kind of genes it annotates,
how many of each kind, whether its protein-coding genes carry InterPro/PFAM domain evidence,
and how many of its genes overlap another gene (same strand / opposite strand), at three
levels: any overlap at all (raw "gene" line span, UTR included), core-span overlap (CDS/exon
extent, UTR excluded), and actual CDS-to-CDS collision.

The raw "gene" span includes any annotated UTR, and UTR-to-UTR overlap between two convergent
genes (pointing at each other, common in compact intergenic regions) is ordinary biology, not
necessarily a sign of a wrong gene call - checked on this project's own A. fumigatus data, 82%
of opposite-strand raw-gene-span overlaps turned out to be UTR-only, with no CDS or even
core-span overlap at all once UTR was excluded. The core-span and CDS-collision counts are
there specifically to see past that UTR noise; the raw gene-span count is kept too since "do
these two gene calls overlap at all" is still useful on its own.

Different sources describe genes differently - one might call a protein-coding gene
"protein_coding_gene" directly (VEuPathDB/FungiDB-style), another plain "gene" with a
gene_biotype=protein_coding attribute (NCBI/GenBank-style), another neither, in which case
this script falls back to looking at the feature's own children (a CDS child means
protein-coding; a tRNA child means a tRNA gene; ...). Every gene this can't confidently
classify goes in the "unclassified" bucket and is also listed in feature_attribute_table.tsv,
so an unfamiliar convention gets caught and checked by hand instead of silently miscounted.

Needs matplotlib (`pip install matplotlib`) for the plots - everything else is standard
library.

Produces, in --out-dir:
    feature_attribute_table.tsv   every distinct (genome, feature type, how it was classified)
                                   combination found and how many genes fell into it - check
                                   this first if a genome's numbers below look wrong
    gene_biotype_counts.tsv       per genome: how many protein-coding (split by domain
                                   evidence present/absent), pseudogene, tRNA, rRNA, other
                                   non-coding, and unclassified genes
    overlap_counts.tsv            per genome: overlapping protein-coding/pseudogene gene
                                   pairs, split by same-strand/opposite-strand, for the raw
                                   gene-span count, the core-span count, and the CDS-collision
                                   count - plus, for the CDS-collision pairs, a
                                   minor/substantial/contained severity breakdown (both strands
                                   combined)
    <genome>_gene_classification.tsv   one row per gene - id, feature type, category, how it
                                   was classified, domain evidence, coordinates
    cds_overlapping_pairs.tsv     one row per CDS-overlapping gene pair, across every genome -
                                   both genes' id, description (from "description=" on the gene
                                   line, or "product=" from a descendant CDS/mRNA as a
                                   fallback), coordinates and CDS length, plus the overlap in bp
                                   and its minor/substantial/contained severity
    01_gene_counts_protein_coding.png  protein-coding genes (split: has domain evidence /
                                   doesn't) + pseudogenes, per genome, totals on top
    02_gene_biotype_breakdown.png      every gene category found, per genome, totals on top
    03_gene_overlap_counts.png         any overlap at all (raw "gene" line span, UTR included),
                                   per genome, split same/opposite strand, totals on top
    04_core_overlap_counts.png         core-span (CDS/exon extent, UTR excluded) overlapping
                                   gene pairs, per genome, split same/opposite strand, totals
                                   on top - a subset of 03's counts, usually much smaller
    05_cds_overlap_counts.png          actual CDS-to-CDS collisions, per genome, split
                                   same/opposite strand, totals on top - a subset of 04's
                                   counts, shown separately since stacking all three together
                                   would make the smaller ones invisible
    06_cds_overlap_severity.png        the same CDS-collision pairs as 05, both strands
                                   combined, split by how much of the smaller gene's own CDS
                                   the overlap covers: minor (<20%), substantial (20-80%),
                                   contained (>=80%, one gene essentially sitting inside the
                                   other) - same rule as step 05's own
                                   plot_final_annotation_qc.py, used later on the post-Liftoff
                                   build, so the two are directly comparable: this plot is the
                                   baseline severity already present before Liftoff runs at
                                   all, so the difference between this and the post-Liftoff
                                   version shows how much of that problem Liftoff introduced
                                   versus how much was already there

Usage:
    python3 analyze_input_gffs.py \\
        --genomes genome_order.txt \\
        --gff-dir gff_cleaned \\
        --out-dir annotation_check
"""

import argparse
import csv
import os
import re
import sys
from urllib.parse import unquote
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
except ImportError:
    sys.exit("ERROR: this script needs matplotlib - install it with 'pip install matplotlib'")

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
GENE_BIOTYPE_RE = re.compile(r"(?:^|;)gene_biotype=([^;]+)")
BIOTYPE_RE = re.compile(r"(?:^|;)biotype=([^;]+)")
EBI_BIOTYPE_RE = re.compile(r"(?:^|;)ebi_biotype=([^;]+)")  # VEuPathDB/FungiDB-style
DESCRIPTION_RE = re.compile(r"(?:^|;)description=([^;]+)")
PRODUCT_RE = re.compile(r"(?:^|;)product=([^;]+)")

# Gene-level feature types whose category is obvious from the type name alone - no need to
# look at attributes or children. "gene" is deliberately not here: it's used by too many
# sources for every biotype, so it always falls through to the attribute/children checks below.
FEATURE_TYPE_TO_CATEGORY = {
    "protein_coding_gene": "protein_coding",
    "pseudogene": "pseudogene",
    "trna_gene": "tRNA",
    "rrna_gene": "rRNA",
    "ncrna_gene": "ncRNA_other",
    "lncrna_gene": "ncRNA_other",
    "snrna_gene": "ncRNA_other",
    "snorna_gene": "ncRNA_other",
    "mirna_gene": "ncRNA_other",
    "scrna_gene": "ncRNA_other",
}

# gene_biotype=/biotype= attribute values (NCBI/GenBank and Ensembl-style conventions).
BIOTYPE_VALUE_TO_CATEGORY = {
    "protein_coding": "protein_coding",
    "mrna": "protein_coding",
    "pseudogene": "pseudogene",
    "trna": "tRNA",
    "rrna": "rRNA",
    "ncrna": "ncRNA_other",
    "lncrna": "ncRNA_other",
    "lnc_rna": "ncRNA_other",
    "snrna": "ncRNA_other",
    "snorna": "ncRNA_other",
    "mirna": "ncRNA_other",
    "misc_rna": "ncRNA_other",
}

# When neither the feature type nor a biotype attribute gives an answer, this is the order
# children types are checked in - CDS first, since a CDS child is definitive regardless of
# whatever else is nested under the same gene.
NCRNA_DESCENDANT_TYPES = {"ncrna", "lnc_rna", "lncrna", "snrna", "snorna", "mirna", "misc_rna"}



def biotype_category(value):
    """Category for a gene_biotype=/biotype= value: exact matches first, then the same idea for the
    many variants other sources use (processed_pseudogene, Mt_tRNA, rRNA_pseudogene, sRNA, ...).
    Anything else (e.g. IG_C_gene) is left unrecognized."""
    v = value.strip().strip('"').lower()
    if v in BIOTYPE_VALUE_TO_CATEGORY:
        return BIOTYPE_VALUE_TO_CATEGORY[v]
    if "pseudogene" in v:
        return "pseudogene"
    if "trna" in v:
        return "tRNA"
    if "rrna" in v:
        return "rRNA"
    if "rna" in v:
        return "ncRNA_other"
    return None


CATEGORIES = ["protein_coding", "pseudogene", "tRNA", "rRNA", "ncRNA_other", "unclassified"]


def read_genomes(value):
    if os.path.isfile(value):
        with open(value) as fh:
            return [l.split()[0] for l in fh if l.strip() and not l.startswith("#")]
    return [g.strip() for g in value.split(",") if g.strip()]


def find_gff(gff_dir, genome):
    for ext in ("gff3", "gff"):
        path = os.path.join(gff_dir, f"{genome}.{ext}")
        if os.path.isfile(path):
            return path
    return None


def has_domain_evidence(attrs):
    return "InterPro:" in attrs or "PFAM:" in attrs or "Pfam:" in attrs


def parse_gff(path):
    """One pass to record every feature; a second walk (via the Parent graph built here)
    to classify each gene-level feature and collect its CDS segments/domain evidence."""
    id_to_feature = {}       # id -> dict(type, start, end, strand, seqname, attrs)
    children_by_parent = defaultdict(list)
    gene_ids = []             # gene-level features, in file order

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqname, _source, ftype, start, end, _score, strand = fields[:7]
            attrs = fields[8]

            id_match = ID_RE.search(attrs)
            feat_id = id_match.group(1).strip() if id_match else None

            if feat_id:
                id_to_feature[feat_id] = {
                    "type": ftype, "start": int(start), "end": int(end),
                    "strand": strand, "seqname": seqname, "attrs": attrs,
                }
                if "gene" in ftype.lower():
                    gene_ids.append(feat_id)

            parent_match = PARENT_RE.search(attrs)
            if parent_match and feat_id:
                for parent_id in parent_match.group(1).split(","):
                    children_by_parent[parent_id.strip()].append(feat_id)

    def walk_descendants(root_id):
        """BFS below one gene: every descendant feature type seen, its CDS segments (start,
        end), whether any CDS carries InterPro/PFAM evidence, and the first "product=" value
        found (NCBI-style genomes put the gene's description there instead of on the gene line
        itself). Falls back to exon segments if the gene has no CDS descendant at all (e.g. a
        pseudogene with no coding sequence annotated) so overlap detection still has something
        to compare."""
        seen_types = set()
        cds_segments, exon_segments = [], []
        domain_hit = False
        product = None
        stack, visited = list(children_by_parent.get(root_id, [])), set()
        while stack:
            cid = stack.pop()
            if cid in visited:
                continue
            visited.add(cid)
            feat = id_to_feature.get(cid)
            if feat is None:
                continue
            ftype_lower = feat["type"].lower()
            seen_types.add(ftype_lower)
            if ftype_lower == "cds":
                cds_segments.append((feat["start"], feat["end"]))
                if has_domain_evidence(feat["attrs"]):
                    domain_hit = True
            elif ftype_lower == "exon":
                exon_segments.append((feat["start"], feat["end"]))
            if product is None:
                m = PRODUCT_RE.search(feat["attrs"])
                if m:
                    product = m.group(1).strip()
            stack.extend(children_by_parent.get(cid, []))
        return seen_types, (cds_segments or exon_segments), domain_hit, product

    genes = []
    attribute_table = defaultdict(int)  # (feature_type, category_source) -> n

    for gene_id in gene_ids:
        feat = id_to_feature[gene_id]
        ftype_lower = feat["type"].lower()
        descendant_types, coding_segments, domain_hit, product = walk_descendants(gene_id)

        desc_match = DESCRIPTION_RE.search(feat["attrs"])
        raw_description = desc_match.group(1).strip() if desc_match else (product or "")
        # GFF3 attribute values percent-encode reserved characters (comma as %2C, etc.) -
        # decode so the description reads as normal text.
        description = unquote(raw_description)

        category, source = None, None
        if ftype_lower in FEATURE_TYPE_TO_CATEGORY:
            category, source = FEATURE_TYPE_TO_CATEGORY[ftype_lower], "feature_type"
            if category == "protein_coding":
                # a protein_coding_gene whose biotype attribute says pseudogene is a pseudogene:
                # this is how the pipeline marks them once gene types are unified (prepare_for_liftoff.py,
                # and the filter step relabelling genes with no valid ORF) - type alone would hide them
                for rx in (GENE_BIOTYPE_RE, BIOTYPE_RE, EBI_BIOTYPE_RE):
                    bm = rx.search(feat["attrs"])
                    if bm and biotype_category(bm.group(1)) == "pseudogene":
                        category, source = "pseudogene", "feature_type + biotype=pseudogene"
                    if bm:
                        break
        else:
            m = GENE_BIOTYPE_RE.search(feat["attrs"])
            attr_name = "gene_biotype"
            if not m:
                m = BIOTYPE_RE.search(feat["attrs"])
                attr_name = "biotype"
            if not m:
                m = EBI_BIOTYPE_RE.search(feat["attrs"])
                attr_name = "ebi_biotype"
            if m:
                value = m.group(1).strip().strip('"').lower()
                category = biotype_category(value)
                source = f"{attr_name}={value}" if category else f"{attr_name}={value} (unrecognized value)"
                category = category or "unclassified"
        if category is None:
            if "cds" in descendant_types:
                category, source = "protein_coding", "descendant:CDS"
            elif "trna" in descendant_types:
                category, source = "tRNA", "descendant:tRNA"
            elif "rrna" in descendant_types:
                category, source = "rRNA", "descendant:rRNA"
            elif descendant_types & NCRNA_DESCENDANT_TYPES:
                category, source = "ncRNA_other", "descendant:ncRNA-like"
            elif "pseudogenic_transcript" in descendant_types or "pseudogenic_exon" in descendant_types:
                category, source = "pseudogene", "descendant:pseudogenic"
            else:
                category, source = "unclassified", "no evidence found"

        # Core span (CDS/exon extent, UTR excluded) used for overlap detection - the raw
        # "gene" line's own start/end includes any annotated UTR, and UTR-to-UTR overlap
        # between convergent genes is common and not itself a sign of a wrong gene call,
        # so it would swamp a gene-span-based overlap count with noise. Falls back to the
        # gene's own start/end only if it has no CDS/exon segments at all (e.g. a
        # pseudogene with no coding sequence annotated at all).
        if coding_segments:
            core_start = min(s for s, _ in coding_segments)
            core_end = max(e for _, e in coding_segments)
        else:
            core_start, core_end = feat["start"], feat["end"]

        attribute_table[(feat["type"], source)] += 1
        genes.append({
            "id": gene_id, "feature_type": feat["type"], "category": category, "source": source,
            "has_domain": domain_hit if category == "protein_coding" else None,
            "description": description,
            "seqname": feat["seqname"], "start": feat["start"], "end": feat["end"],
            "core_start": core_start, "core_end": core_end,
            "strand": feat["strand"], "cds_segments": coding_segments,
        })

    return genes, attribute_table


def segments_overlap_bp(segs1, segs2):
    total = 0
    for s1, e1 in segs1:
        for s2, e2 in segs2:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi >= lo:
                total += hi - lo + 1
    return total


def cds_bp(gene):
    return sum(e - s + 1 for s, e in gene["cds_segments"])


def overlap_severity(overlap_bp, gene_a, gene_b):
    """Same rule as this pipeline's later plot_final_annotation_qc.py (step 05's own QC on
    the post-Liftoff build) - severity is the overlap as a fraction of the SMALLER gene's own
    CDS length, not the bigger one and not the two combined. Sized against the smaller gene on
    purpose: a 438bp overlap is only 22% of a 1986bp gene, but 100% of a 438bp one - sizing
    against the bigger gene would call that pair "minor" even though the smaller gene is
    entirely swallowed by the overlap, which is exactly the case most worth flagging."""
    smaller_bp = min(cds_bp(gene_a), cds_bp(gene_b))
    if smaller_bp == 0:
        return None
    pct = overlap_bp / smaller_bp
    if pct >= 0.8:
        return "contained"
    elif pct >= 0.2:
        return "substantial"
    else:
        return "minor"


def find_overlaps(genes, start_key, end_key):
    """Sweep-line scan (same approach as check_overlapping_genes.py) restricted to
    protein-coding and pseudogene genes - overlap between two non-coding RNA genes, or between
    a non-coding gene and a coding one, isn't the same-locus-consistency problem this pipeline
    cares about here. `start_key`/`end_key` pick which span to sweep on:
    "start"/"end" for the raw "gene" line's own coordinates (any overlap at all, UTR
    included - "does the gene, from its very first to very last base, overlap another
    gene's"), or "core_start"/"core_end" for the CDS/exon extent (UTR excluded). Either way,
    `cds_overlap` on each returned pair is the actual CDS-to-CDS collision - a pair can have a
    "gene" or "core" overlap without a CDS overlap, but not the other way around."""
    relevant = [g for g in genes if g["category"] in ("protein_coding", "pseudogene")]
    by_seq = defaultdict(list)
    for gene in relevant:
        by_seq[gene["seqname"]].append(gene)

    overlaps = []
    for seq_genes in by_seq.values():
        seq_genes.sort(key=lambda g: g[start_key])
        active = []
        for gene in seq_genes:
            active = [a for a in active if a[end_key] >= gene[start_key]]
            for other in active:
                ov_bp = segments_overlap_bp(other["cds_segments"], gene["cds_segments"])
                overlaps.append({
                    "gene_a": other, "gene_b": gene, "overlap_bp": ov_bp,
                    "same_strand": other["strand"] == gene["strand"],
                    "cds_overlap": ov_bp > 0,
                    "severity": overlap_severity(ov_bp, other, gene) if ov_bp > 0 else None,
                })
            active.append(gene)
    return overlaps


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genomes", required=True, help="Comma-separated genome names or a file with one per line (e.g. genome_order.txt)")
    parser.add_argument("--gff-dir", required=True, help="Directory with each genome's own <genome>.gff3 or <genome>.gff")
    parser.add_argument("--out-dir", default="./annotation_check", help="Where to write the tables and plots")
    args = parser.parse_args()

    genomes = read_genomes(args.genomes)
    os.makedirs(args.out_dir, exist_ok=True)

    biotype_rows, overlap_rows, pair_rows = [], [], []
    combined_attribute_table = defaultdict(int)  # (genome, feature_type, source) -> n

    for genome in genomes:
        gff_path = find_gff(args.gff_dir, genome)
        if gff_path is None:
            print(f"WARNING: no GFF found for {genome} in {args.gff_dir} (tried .gff3, .gff) - skipping", file=sys.stderr)
            continue

        print(f"[{genome}] loading {gff_path}...", file=sys.stderr)
        genes, attribute_table = parse_gff(gff_path)
        if not genes:
            print(f"WARNING: [{genome}] 0 gene-level features found - this script only "
                  f"recognizes a feature type as gene-level if its name contains \"gene\" "
                  f"(gene, protein_coding_gene, pseudogene, tRNA_gene, ...); if this genome's "
                  f"GFF uses something else, this genome will be reported as all-zero below",
                  file=sys.stderr)
        for (ftype, source), n in attribute_table.items():
            combined_attribute_table[(genome, ftype, source)] += n

        with open(os.path.join(args.out_dir, f"{genome}_gene_classification.tsv"), "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["gene_id", "feature_type", "category", "classified_by", "has_domain", "seqname", "start", "end", "strand"])
            for g in genes:
                writer.writerow([g["id"], g["feature_type"], g["category"], g["source"],
                                  "" if g["has_domain"] is None else g["has_domain"],
                                  g["seqname"], g["start"], g["end"], g["strand"]])

        counts = {cat: 0 for cat in CATEGORIES}
        pc_with_domain = pc_without_domain = 0
        for g in genes:
            counts[g["category"]] += 1
            if g["category"] == "protein_coding":
                if g["has_domain"]:
                    pc_with_domain += 1
                else:
                    pc_without_domain += 1

        n_unclassified = counts["unclassified"]
        if n_unclassified:
            print(f"WARNING: [{genome}] {n_unclassified} gene(s) could not be classified - "
                  f"see feature_attribute_table.tsv", file=sys.stderr)

        biotype_rows.append({
            "genome": genome,
            "protein_coding": counts["protein_coding"],
            "protein_coding_with_domain": pc_with_domain,
            "protein_coding_without_domain": pc_without_domain,
            "pseudogene": counts["pseudogene"],
            "tRNA": counts["tRNA"], "rRNA": counts["rRNA"],
            "ncRNA_other": counts["ncRNA_other"], "unclassified": counts["unclassified"],
            "total_genes": len(genes),
        })

        gene_overlaps = find_overlaps(genes, "start", "end")
        same_gene = sum(1 for ov in gene_overlaps if ov["same_strand"])
        opp_gene = sum(1 for ov in gene_overlaps if not ov["same_strand"])

        core_overlaps = find_overlaps(genes, "core_start", "core_end")
        same_core = sum(1 for ov in core_overlaps if ov["same_strand"])
        same_cds = sum(1 for ov in core_overlaps if ov["same_strand"] and ov["cds_overlap"])
        opp_core = sum(1 for ov in core_overlaps if not ov["same_strand"])
        opp_cds = sum(1 for ov in core_overlaps if not ov["same_strand"] and ov["cds_overlap"])

        severity_counts = {"minor": 0, "substantial": 0, "contained": 0}
        for ov in core_overlaps:
            if ov["severity"] is not None:
                severity_counts[ov["severity"]] += 1
                a, b = ov["gene_a"], ov["gene_b"]
                pair_rows.append({
                    "genome": genome, "seqname": a["seqname"],
                    "same_strand": ov["same_strand"],
                    "gene1_id": a["id"], "gene1_description": a["description"],
                    "gene1_start": a["start"], "gene1_end": a["end"], "gene1_cds_bp": cds_bp(a),
                    "gene1_has_domain": "" if a["has_domain"] is None else a["has_domain"],
                    "gene2_id": b["id"], "gene2_description": b["description"],
                    "gene2_start": b["start"], "gene2_end": b["end"], "gene2_cds_bp": cds_bp(b),
                    "gene2_has_domain": "" if b["has_domain"] is None else b["has_domain"],
                    "overlap_bp": ov["overlap_bp"], "severity": ov["severity"],
                })

        overlap_rows.append({
            "genome": genome,
            "same_strand_gene_overlap": same_gene, "opposite_strand_gene_overlap": opp_gene,
            "same_strand_core_overlap": same_core, "opposite_strand_core_overlap": opp_core,
            "same_strand_cds_overlap": same_cds, "opposite_strand_cds_overlap": opp_cds,
            "cds_overlap_minor": severity_counts["minor"],
            "cds_overlap_substantial": severity_counts["substantial"],
            "cds_overlap_contained": severity_counts["contained"],
        })

    if not biotype_rows:
        sys.exit("ERROR: no genome had a GFF to analyze - nothing done")

    with open(os.path.join(args.out_dir, "gene_biotype_counts.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(biotype_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(biotype_rows)

    with open(os.path.join(args.out_dir, "overlap_counts.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(overlap_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(overlap_rows)

    with open(os.path.join(args.out_dir, "cds_overlapping_pairs.tsv"), "w", newline="") as fh:
        fieldnames = ["genome", "seqname", "same_strand", "severity", "overlap_bp",
                      "gene1_id", "gene1_description", "gene1_start", "gene1_end", "gene1_cds_bp", "gene1_has_domain",
                      "gene2_id", "gene2_description", "gene2_start", "gene2_end", "gene2_cds_bp", "gene2_has_domain"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_rows)

    with open(os.path.join(args.out_dir, "feature_attribute_table.tsv"), "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["genome", "feature_type", "classified_by", "n"])
        for (genome, ftype, source), n in sorted(combined_attribute_table.items()):
            writer.writerow([genome, ftype, source, n])

    analyzed_genomes = [r["genome"] for r in biotype_rows]
    figsize = (max(8, len(analyzed_genomes) * 0.6), 5)

    def add_totals(ax, x_positions, totals):
        for xpos, total in zip(x_positions, totals):
            ax.text(xpos, total, str(total), ha="center", va="bottom", fontsize=8)
        # Counts are always whole numbers - force integer-only y-axis ticks (matplotlib's
        # default locator can otherwise pick fractional ticks like 2.5, 7.5 for small ranges).
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # Headroom so the total label above the tallest bar doesn't sit on/above the plot's
        # own top border - grows the ylim, never shrinks it, so this is safe to call more than
        # once per plot (once per same/opposite-strand series).
        current_top = ax.get_ylim()[1]
        needed_top = max(totals, default=0) * 1.15
        ax.set_ylim(top=max(current_top, needed_top, 1))

    # --- Plot 1: protein-coding genes (split by domain evidence) + pseudogenes ---
    fig, ax = plt.subplots(figsize=figsize)
    x = list(range(len(analyzed_genomes)))
    with_domain = [r["protein_coding_with_domain"] for r in biotype_rows]
    without_domain = [r["protein_coding_without_domain"] for r in biotype_rows]
    pseudo = [r["pseudogene"] for r in biotype_rows]
    ax.bar(x, with_domain, label="protein-coding, has domain evidence", color="#2F7A55")
    ax.bar(x, without_domain, bottom=with_domain, label="protein-coding, no domain evidence", color="#B0B0B0")
    bottoms2 = [w + n for w, n in zip(with_domain, without_domain)]
    ax.bar(x, pseudo, bottom=bottoms2, label="pseudogene", color="#C0392B")
    totals1 = [w + n + p for w, n, p in zip(with_domain, without_domain, pseudo)]
    add_totals(ax, x, totals1)
    ax.set_xticks(x)
    ax.set_xticklabels(analyzed_genomes, rotation=45, ha="right")
    ax.set_ylabel("Genes")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "01_gene_counts_protein_coding.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 2: full gene biotype breakdown ---
    fig, ax = plt.subplots(figsize=figsize)
    breakdown_colors = {
        "protein_coding": "#3B7DD8", "pseudogene": "#C0392B", "tRNA": "#F2B035",
        "rRNA": "#B379E0", "ncRNA_other": "#1ABC9C", "unclassified": "#7A5195",
    }
    bottoms = [0] * len(analyzed_genomes)
    for cat in CATEGORIES:
        values = [r[cat] for r in biotype_rows]
        ax.bar(x, values, bottom=bottoms, label=cat, color=breakdown_colors[cat])
        bottoms = [b + v for b, v in zip(bottoms, values)]
    add_totals(ax, x, bottoms)
    ax.set_xticks(x)
    ax.set_xticklabels(analyzed_genomes, rotation=45, ha="right")
    ax.set_ylabel("Genes")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "02_gene_biotype_breakdown.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 3: raw gene-span overlap counts (any overlap at all, UTR included) ---
    fig, ax = plt.subplots(figsize=(max(10, len(analyzed_genomes) * 1.0), 5))
    width = 0.35
    x_same = [i - width / 2 for i in x]
    x_opp = [i + width / 2 for i in x]
    same_gene = [r["same_strand_gene_overlap"] for r in overlap_rows]
    opp_gene = [r["opposite_strand_gene_overlap"] for r in overlap_rows]
    ax.bar(x_same, same_gene, width, label="same strand", color="#B0B0B0")
    ax.bar(x_opp, opp_gene, width, label="opposite strand", color="#619CFF")
    add_totals(ax, x_same, same_gene)
    add_totals(ax, x_opp, opp_gene)
    ax.set_xticks(x)
    ax.set_xticklabels(analyzed_genomes, rotation=45, ha="right")
    ax.set_ylabel("Overlapping gene pairs (raw gene span - UTR included)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "03_gene_overlap_counts.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 4: core-span overlap counts (CDS/exon extent, UTR excluded) ---
    fig, ax = plt.subplots(figsize=(max(10, len(analyzed_genomes) * 1.0), 5))
    same_core = [r["same_strand_core_overlap"] for r in overlap_rows]
    opp_core = [r["opposite_strand_core_overlap"] for r in overlap_rows]
    ax.bar(x_same, same_core, width, label="same strand", color="#B0B0B0")
    ax.bar(x_opp, opp_core, width, label="opposite strand", color="#619CFF")
    add_totals(ax, x_same, same_core)
    add_totals(ax, x_opp, opp_core)
    ax.set_xticks(x)
    ax.set_xticklabels(analyzed_genomes, rotation=45, ha="right")
    ax.set_ylabel("Overlapping gene pairs (core span - UTR excluded)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "04_core_overlap_counts.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 5: actual CDS-to-CDS collisions ---
    fig, ax = plt.subplots(figsize=(max(10, len(analyzed_genomes) * 1.0), 5))
    same_cds = [r["same_strand_cds_overlap"] for r in overlap_rows]
    opp_cds = [r["opposite_strand_cds_overlap"] for r in overlap_rows]
    ax.bar(x_same, same_cds, width, label="same strand", color="#C0392B")
    ax.bar(x_opp, opp_cds, width, label="opposite strand", color="#F2B035")
    add_totals(ax, x_same, same_cds)
    add_totals(ax, x_opp, opp_cds)
    ax.set_xticks(x)
    ax.set_xticklabels(analyzed_genomes, rotation=45, ha="right")
    ax.set_ylabel("Overlapping gene pairs (CDS-to-CDS collision)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "05_cds_overlap_counts.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 6: CDS-overlap severity (same rule as step 05's plot_final_annotation_qc.py) ---
    fig, ax = plt.subplots(figsize=figsize)
    sev_colors = {"minor": "#B0B0B0", "substantial": "#F2B035", "contained": "#C0392B"}
    bottoms = [0] * len(analyzed_genomes)
    for sev in ("minor", "substantial", "contained"):
        values = [r[f"cds_overlap_{sev}"] for r in overlap_rows]
        ax.bar(x, values, bottom=bottoms, label=sev, color=sev_colors[sev])
        bottoms = [b + v for b, v in zip(bottoms, values)]
    add_totals(ax, x, bottoms)
    ax.set_xticks(x)
    ax.set_xticklabels(analyzed_genomes, rotation=45, ha="right")
    ax.set_ylabel("CDS-overlapping pairs (both strands)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "06_cds_overlap_severity.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("\nGene biotype counts:", file=sys.stderr)
    header = ["genome"] + CATEGORIES + ["total_genes"]
    print("\t".join(header), file=sys.stderr)
    for r in biotype_rows:
        print("\t".join(str(r[c]) for c in header), file=sys.stderr)

    print(f"\nDone - tables + 6 plots written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
