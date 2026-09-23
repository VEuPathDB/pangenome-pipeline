#!/usr/bin/env python3
"""
QC plots for this pipeline's raw output - <genome>_final_annotation.gff. Same core plots
this project has been using all along (gene counts, valid-ORF rate, extra copies, overlap
severity, source composition), just as a standalone Python/matplotlib script instead of an
R Markdown report - run this as a first look at the raw build's quality. No filtering step
exists yet in this folder - this script only reads and reports, it never changes any GFF.

Self-contained - the GFF-parsing and CDS-overlap-detection logic below is its own copy, not
shared with anything else in this folder.

Needs matplotlib (`pip install matplotlib`) - everything else is standard library.

Produces, in --out-dir:
    01_gene_counts.png       protein-coding genes: own original annotation vs. this pipeline's output
    02_valid_orf_rate.png    % of genes with a valid ORF, per strain
    03_extra_copies.png      extra-copy genes, split by valid ORF, per strain
    04_overlap_severity.png  CDS-overlap severity (minor/substantial/contained), per strain
    05_source_composition.png  where the genes came from (main/original/good_*/rescued_*), per strain
    06_overlap_strand.png      the same pairs as 04, split same-strand vs opposite-strand instead
                              of by severity - same-strand pairs are the ones "that shouldn't
                              normally happen"
    07_overlap_domain_support.png  the same pairs as 04, split by whether they have InterPro/PFAM
                              domain evidence on both genes, one, or neither
    08_minor_overlap_position.png  minor-severity pairs only (domain evidence not considered here),
                              split by where the overlap sits in the bigger gene's own exon
                              structure - flank (an outermost exon - likely just boundary
                              imprecision, both genes real) vs inside (an internal exon, meaning
                              the smaller gene mostly sits in an intron - likely spurious) vs
                              single-exon (bigger gene has no intron to be "inside" of)
    summary.tsv              the numbers behind every plot above, one strain per row
    extra_copy_genes.tsv     every extra-copy gene, one row per strain: id, description, location,
                              source, valid ORF, Liftoff coverage/identity - for looking up what
                              these genes actually are, not just how many
    cds_overlapping_pairs.tsv  one row per CDS-overlapping gene pair, across every strain - same
                              columns as this project's pre-Liftoff annotation check
                              (analyze_input_gffs.py), plus each gene's own source layer
    overlap_domain_summary.tsv  per strain: how many of the overlapping pairs above have InterPro/
                              PFAM domain evidence on both genes, on one, or on neither, and how
                              many DISTINCT genes are involved in at least one overlap (a gene in
                              3 pairs only counts once)

Usage:
    python3 plot_final_annotation_qc.py \\
        --strains "AF293,strain2,strain3" \\
        --final-dir ./expanded_annotation \\
        --genome-dir /path/to/each_genomes_own_annotation \\
        --genome-suffix "_genomic.gff" \\
        --out-dir ./final_annotation_qc_plots
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from urllib.parse import unquote

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("ERROR: this script needs matplotlib - install it with 'pip install matplotlib'")

GENE_BIOTYPE_RE = re.compile(r"(?:^|;)gene_biotype=([^;]+)")
ID_RE = re.compile(r"ID=([^;]+)")
PARENT_RE = re.compile(r"Parent=([^;]+)")
DESCRIPTION_RE = re.compile(r"(?:^|;)description=([^;]+)")
PRODUCT_RE = re.compile(r"(?:^|;)product=([^;]+)")
COPY_NUM_ID_RE = re.compile(r"(?:^|;)copy_num_ID=([^;]+)")
COVERAGE_RE = re.compile(r"(?:^|;)coverage=([^;]+)")
SEQUENCE_ID_RE = re.compile(r"(?:^|;)sequence_ID=([^;]+)")


def _int_attr(attrs, key):
    m = re.search(rf"{key}=(\d+)", attrs)
    return int(m.group(1)) if m else None


def _valid_orf(attrs):
    if "valid_ORF=True" in attrs:
        return True
    if "valid_ORF=False" in attrs:
        return False
    m = re.search(r"valid_ORFs=(\d+)", attrs)
    if m:
        return int(m.group(1)) >= 1
    return None


def parse_gff(path):
    """Gene-level records, mRNA->gene links, and per-mRNA CDS segments/domain flags -
    same logic this project's filter_final_annotation.py uses, kept here as its own copy
    so this QC script has no dependency on that (or any other) file."""
    genes = {}       # gene_id -> dict of fields/attrs
    mrna_to_gene = {}
    cds_by_mrna = defaultdict(list)   # mrna_id -> [(start, end)]
    mrna_has_domain = defaultdict(bool)

    with open(path) as fh:
        raw_lines = fh.readlines()

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqname, source, feature, start, end, score, strand = fields[:7]
        attrs = fields[8]
        m = ID_RE.search(attrs)
        feat_id = m.group(1) if m else None

        if feature == "protein_coding_gene" and feat_id:
            desc_m = DESCRIPTION_RE.search(attrs) or PRODUCT_RE.search(attrs)
            cov_m, seqid_m, copy_m = COVERAGE_RE.search(attrs), SEQUENCE_ID_RE.search(attrs), COPY_NUM_ID_RE.search(attrs)
            genes[feat_id] = {
                "id": feat_id, "seqname": seqname, "source": source,
                "start": int(start), "end": int(end), "strand": strand,
                "extra_copy_number": _int_attr(attrs, "extra_copy_number"),
                "valid_ORF": _valid_orf(attrs),
                "description": unquote(desc_m.group(1)) if desc_m else "",
                "coverage": cov_m.group(1) if cov_m else "",
                "sequence_ID": seqid_m.group(1) if seqid_m else "",
                "copy_num_ID": copy_m.group(1) if copy_m else "",
            }
        elif feature == "mRNA" and feat_id:
            parent = PARENT_RE.search(attrs)
            if parent:
                mrna_to_gene[feat_id] = parent.group(1)

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        feature, start, end, attrs = fields[2], fields[3], fields[4], fields[8]
        if feature != "CDS":
            continue
        parent = PARENT_RE.search(attrs)
        if not parent:
            continue
        mrna_id = parent.group(1)
        if mrna_to_gene.get(mrna_id) not in genes:
            continue
        cds_by_mrna[mrna_id].append((int(start), int(end)))
        if "InterPro:" in attrs or "PFAM:" in attrs or "Pfam:" in attrs:
            mrna_has_domain[mrna_id] = True

    # Longest-isoform CDS segments and domain flag, per gene.
    for gene_id, gene in genes.items():
        gene_mrnas = [m for m, g in mrna_to_gene.items() if g == gene_id]
        best_mrna, best_bp = None, -1
        for mrna_id in gene_mrnas:
            segs = cds_by_mrna.get(mrna_id, [])
            bp = sum(e - s + 1 for s, e in segs)
            if bp > best_bp:
                best_bp, best_mrna = bp, mrna_id
        # sorted left-to-right by genomic position, not transcript/exon order - "outermost" here
        # always means the segment with the lowest start or the highest end, regardless of strand.
        gene["cds_segments"] = sorted(cds_by_mrna.get(best_mrna, []))
        gene["cds_bp"] = best_bp if best_bp > 0 else 0
        gene["has_domain"] = any(mrna_has_domain.get(m, False) for m in gene_mrnas)

    return genes


def segments_overlap_bp(segs1, segs2):
    total = 0
    for s1, e1 in segs1:
        for s2, e2 in segs2:
            lo, hi = max(s1, s2), min(e1, e2)
            if hi >= lo:
                total += hi - lo + 1
    return total


def classify_pair_case(g1, g2):
    """Which of the 6 mutually-exclusive overlap-resolution cases this pair falls into (see the
    "Refined overlap-resolution rule" in the general README, step 7) - computed the same way
    there: is either gene `main`, and does each carry domain evidence."""
    m1, m2 = g1["source"] == "main", g2["source"] == "main"
    d1, d2 = g1["has_domain"], g2["has_domain"]
    if m1 and m2:
        return "1_main_x_main"
    if m1 or m2:
        other_has_domain = d2 if m1 else d1
        return "2_main_x_other_with_domain" if other_has_domain else "3_main_x_other_no_domain"
    if d1 != d2:
        return "4_other_x_other_one_domain"
    return "5_other_x_other_both_domain" if d1 else "6_other_x_other_no_domain"


def classify_minor_position(bigger_gene, ov_bp, smaller_cds_segments):
    """Only meaningful for a "minor"-severity pair: does the overlap sit at one of the BIGGER
    gene's own two outermost CDS segments ("flank" - the overlap is near where the bigger gene's
    own model could legitimately be a bit off, or where a genuinely separate adjacent gene could
    sit just outside its real boundary - most likely just gene-model boundary imprecision, both
    genes real), or on one of its INTERNAL segments ("inside" - the smaller gene's CDS reaches an
    exon buried in the middle of the bigger gene, meaning most of the smaller gene's own body
    sits inside one of the bigger gene's introns - much harder to explain as boundary noise, more
    likely a spurious gene call for the smaller one). "single_exon" when the bigger gene has only
    one CDS segment - the flank/inside distinction needs at least one intron to mean anything."""
    segs = bigger_gene["cds_segments"]
    if len(segs) <= 1:
        return "single_exon"
    outer = {segs[0], segs[-1]}
    overlapped = {s for s in segs for ss, se in smaller_cds_segments
                  if max(s[0], ss) <= min(s[1], se)}
    if not overlapped:
        return None  # shouldn't happen if ov_bp > 0, but don't crash on it
    return "flank" if overlapped <= outer else "inside"


def find_cds_overlaps(genes):
    """Sweep-line scan for CDS-to-CDS overlapping gene pairs, returning the overlap length
    in bp and whether the pair is same-strand for each pair - same approach as
    check_overlapping_genes.py and analyze_input_gffs.py's pre-Liftoff version of this same
    check (identical segments_overlap_bp, identical "smaller gene's CDS" denominator, identical
    20%/80% thresholds - the two are meant to be directly comparable, see severity_counts)."""
    by_seq = defaultdict(list)
    for gene in genes.values():
        by_seq[gene["seqname"]].append(gene)

    pairs = []
    for seqname, seq_genes in by_seq.items():
        seq_genes.sort(key=lambda g: g["start"])
        active = []
        for gene in seq_genes:
            active = [a for a in active if a["end"] >= gene["start"]]
            for other in active:
                ov_bp = segments_overlap_bp(other["cds_segments"], gene["cds_segments"])
                if ov_bp > 0:
                    pairs.append((other["id"], gene["id"], ov_bp, other["strand"] == gene["strand"]))
            active.append(gene)
    return pairs


def read_strains(value):
    if os.path.isfile(value):
        with open(value) as fh:
            return [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    return [s.strip() for s in value.split(",") if s.strip()]


def count_original_protein_coding(path):
    """Works for either naming convention: FungiDB-style 'protein_coding_gene' as its own
    feature type, or NCBI-style plain 'gene' with a gene_biotype=protein_coding attribute."""
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature, attrs = fields[2], fields[8]
            if feature == "protein_coding_gene":
                n += 1
            elif feature == "gene":
                m = GENE_BIOTYPE_RE.search(attrs)
                if m and m.group(1) == "protein_coding":
                    n += 1
    return n


def classify_severity(ov_bp, cds_bp_a, cds_bp_b):
    """Same rule as analyze_input_gffs.py's overlap_severity(): the overlap as a fraction of the
    SMALLER gene's own CDS length, not the bigger one and not the two combined. Returns None if
    neither gene has a usable CDS length (nothing to size the fraction against)."""
    smaller_bp = min(cds_bp_a, cds_bp_b)
    if smaller_bp == 0:
        return None
    pct = ov_bp / smaller_bp
    if pct >= 0.8:
        return "contained"
    elif pct >= 0.2:
        return "substantial"
    return "minor"


def severity_counts(genes):
    """Same classification as filter_final_annotation.py's overlap-resolution step. Returns the
    minor/substantial/contained counts (both strands combined, same as this project's pre-Liftoff
    baseline in analyze_input_gffs.py) plus a same-strand/opposite-strand split of the raw pair
    count - same-strand CDS collisions are the ones that "shouldn't normally happen" (see
    check_overlapping_genes.py); most of what the all-against-all Liftoff introduces is
    opposite-strand, because AGAT's own complement step only ever rejects a gene for colliding
    with a SAME-strand CDS already placed - an opposite-strand collision, or two genes added in
    the same round that collide with each other, are never blocked at that stage. That is exactly
    what step 7's filter resolves, so a high pair count here, before filtering, is expected."""
    pairs = find_cds_overlaps(genes)
    counts = {"minor": 0, "substantial": 0, "contained": 0}
    strand = {"same_strand": 0, "opposite_strand": 0}
    for gene1_id, gene2_id, ov_bp, same_strand in pairs:
        g1, g2 = genes[gene1_id], genes[gene2_id]
        severity = classify_severity(ov_bp, g1["cds_bp"], g2["cds_bp"])
        if severity is None:
            continue
        strand["same_strand" if same_strand else "opposite_strand"] += 1
        counts[severity] += 1
    counts.update(strand)
    return counts


def source_counts(genes, reference_strain):
    """reference_strain is the actual genome name (e.g. "AF293"), not the "main" tag - needed
    to split "rescued_<reference>" out from every other strain's "rescued_<x>"."""
    ref_rescued_tag = f"rescued_{reference_strain}"
    counts = {"main": 0, "original": 0, "reference_rescued": 0, "good_others": 0, "rescued_others": 0}
    for g in genes.values():
        src = g["source"]
        if src == "main":
            counts["main"] += 1
        elif src == ref_rescued_tag:
            counts["reference_rescued"] += 1
        elif src.startswith("good_"):
            counts["good_others"] += 1
        elif src.startswith("rescued_"):
            counts["rescued_others"] += 1
        else:
            counts["original"] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--strains", required=True, help="Comma-separated genome names or a file with one per line - first one must be the reference")
    parser.add_argument("--final-dir", required=True, help="Directory with <strain>_final_annotation.gff (this pipeline's raw output, before filtering)")
    parser.add_argument("--genome-dir", required=True, help="Directory with each strain's own original GFF")
    parser.add_argument("--genome-suffix", required=True, help="Filename suffix for each strain's own original GFF, e.g. '_genomic.gff' (no default - varies by data source)")
    parser.add_argument("--out-dir", default="./final_annotation_qc_plots", help="Where to write the plots and summary.tsv")
    args = parser.parse_args()

    strains = read_strains(args.strains)
    os.makedirs(args.out_dir, exist_ok=True)

    rows = []
    severity_all, source_all, minor_position_all = {}, {}, {}
    extra_copy_rows = []
    pair_rows = []
    domain_summary_rows = []

    for s in strains:
        final_path = os.path.join(args.final_dir, f"{s}_final_annotation.gff")
        genome_path = os.path.join(args.genome_dir, f"{s}{args.genome_suffix}")
        for p in (final_path, genome_path):
            if not os.path.isfile(p):
                sys.exit(f"ERROR: missing {p}")

        print(f"[{s}] loading...", file=sys.stderr)
        genes = parse_gff(final_path)
        n_original = count_original_protein_coding(genome_path)
        n_final = len(genes)
        n_valid_orf = sum(1 for g in genes.values() if g["valid_ORF"])
        n_extra_copy = sum(1 for g in genes.values() if g["extra_copy_number"])
        n_extra_copy_valid = sum(1 for g in genes.values() if g["extra_copy_number"] and g["valid_ORF"])
        n_extra_copy_invalid = n_extra_copy - n_extra_copy_valid

        sev = severity_counts(genes)
        severity_all[s] = sev
        source_all[s] = source_counts(genes, strains[0])

        # Pair-level detail (mirrors this project's pre-Liftoff cds_overlapping_pairs.tsv) plus
        # per-strain domain-evidence and distinct-gene summaries for those same pairs.
        overlap_pairs = find_cds_overlaps(genes)
        domain_counts = {"both": 0, "one": 0, "none": 0}
        minor_position_counts = {"flank": 0, "inside": 0, "single_exon": 0}
        genes_in_overlap = set()
        for gene1_id, gene2_id, ov_bp, same_strand in overlap_pairs:
            g1, g2 = genes[gene1_id], genes[gene2_id]
            severity = classify_severity(ov_bp, g1["cds_bp"], g2["cds_bp"])
            if severity is None:
                continue
            genes_in_overlap.update((gene1_id, gene2_id))
            n_with_domain = int(g1["has_domain"]) + int(g2["has_domain"])
            domain_counts["both" if n_with_domain == 2 else "one" if n_with_domain == 1 else "none"] += 1
            # Only meaningful for "minor" severity - see classify_minor_position()'s own docstring.
            # Sized against the BIGGER gene on purpose here (opposite of the severity % denominator,
            # which uses the smaller one) - what matters for this classification is whether the
            # overlap sits at the edge or in the middle of the LARGER gene's own exon structure.
            if severity == "minor":
                bigger, smaller = (g1, g2) if g1["cds_bp"] >= g2["cds_bp"] else (g2, g1)
                minor_position = classify_minor_position(bigger, ov_bp, smaller["cds_segments"])
                if minor_position:
                    minor_position_counts[minor_position] += 1
            else:
                minor_position = ""
            mentions_cytochrome = "cytochrome" in g1["description"].lower() or "cytochrome" in g2["description"].lower()
            pair_rows.append({
                "genome": s, "seqname": g1["seqname"], "same_strand": same_strand, "severity": severity,
                "overlap_bp": ov_bp, "minor_position": minor_position, "mentions_cytochrome": mentions_cytochrome,
                "gene1_id": g1["id"], "gene1_description": g1["description"],
                "gene1_start": g1["start"], "gene1_end": g1["end"], "gene1_cds_bp": g1["cds_bp"],
                "gene1_has_domain": g1["has_domain"],
                "gene1_valid_ORF": "NA" if g1["valid_ORF"] is None else g1["valid_ORF"],
                "gene1_source": g1["source"],
                "gene2_id": g2["id"], "gene2_description": g2["description"],
                "gene2_start": g2["start"], "gene2_end": g2["end"], "gene2_cds_bp": g2["cds_bp"],
                "gene2_has_domain": g2["has_domain"],
                "gene2_valid_ORF": "NA" if g2["valid_ORF"] is None else g2["valid_ORF"],
                "gene2_source": g2["source"],
            })
        minor_position_all[s] = minor_position_counts
        n_pairs = sum(domain_counts.values())
        domain_summary_rows.append({
            "strain": s, "n_pairs": n_pairs, "n_genes_in_overlap": len(genes_in_overlap),
            "n_both_domain": domain_counts["both"],
            "pct_both_domain": round(100 * domain_counts["both"] / n_pairs, 1) if n_pairs else 0,
            "n_one_domain": domain_counts["one"],
            "pct_one_domain": round(100 * domain_counts["one"] / n_pairs, 1) if n_pairs else 0,
            "n_no_domain": domain_counts["none"],
            "pct_no_domain": round(100 * domain_counts["none"] / n_pairs, 1) if n_pairs else 0,
        })

        for g in genes.values():
            if not g["extra_copy_number"]:
                continue
            extra_copy_rows.append({
                "strain": s, "gene_id": g["id"], "description": g["description"],
                "seqname": g["seqname"], "start": g["start"], "end": g["end"], "strand": g["strand"],
                "source": g["source"], "extra_copy_number": g["extra_copy_number"],
                "copy_num_ID": g["copy_num_ID"],
                "valid_ORF": "NA" if g["valid_ORF"] is None else g["valid_ORF"],
                "liftoff_coverage": g["coverage"], "liftoff_sequence_ID": g["sequence_ID"],
            })

        rows.append({
            "strain": s, "n_original": n_original, "n_final": n_final,
            "n_valid_orf": n_valid_orf,
            "pct_valid_orf": round(100 * n_valid_orf / n_final, 1) if n_final else 0,
            "n_extra_copy": n_extra_copy,
            "n_extra_copy_valid_orf": n_extra_copy_valid, "n_extra_copy_invalid_orf": n_extra_copy_invalid,
            "n_overlap_pairs": sev["minor"] + sev["substantial"] + sev["contained"],
            "n_minor": sev["minor"], "n_substantial": sev["substantial"], "n_contained": sev["contained"],
            "n_same_strand_overlap": sev["same_strand"], "n_opposite_strand_overlap": sev["opposite_strand"],
            "n_genes_in_overlap": len(genes_in_overlap),
        })

    with open(os.path.join(args.out_dir, "summary.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with open(os.path.join(args.out_dir, "cds_overlapping_pairs.tsv"), "w", newline="") as fh:
        fieldnames = ["genome", "seqname", "same_strand", "severity", "overlap_bp",
                      "minor_position", "mentions_cytochrome",
                      "gene1_id", "gene1_description", "gene1_start", "gene1_end", "gene1_cds_bp",
                      "gene1_has_domain", "gene1_valid_ORF", "gene1_source",
                      "gene2_id", "gene2_description", "gene2_start", "gene2_end", "gene2_cds_bp",
                      "gene2_has_domain", "gene2_valid_ORF", "gene2_source"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_rows)

    with open(os.path.join(args.out_dir, "overlap_domain_summary.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(domain_summary_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(domain_summary_rows)

    # Sorted by strain, then most-duplicated gene first (extra_copy_number 1 = the 2nd copy overall,
    # so the highest number here is the most-repeated locus in that genome) - the busiest loci float
    # to the top of each strain's block for a quick look at what is actually being over-called.
    extra_copy_rows.sort(key=lambda r: (r["strain"], -r["extra_copy_number"]))
    with open(os.path.join(args.out_dir, "extra_copy_genes.tsv"), "w", newline="") as fh:
        fieldnames = ["strain", "gene_id", "description", "seqname", "start", "end", "strand",
                      "source", "extra_copy_number", "copy_num_ID", "valid_ORF",
                      "liftoff_coverage", "liftoff_sequence_ID"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(extra_copy_rows)

    figsize = (max(8, len(strains) * 0.6), 5)
    x = range(len(strains))

    # --- Plot 1: gene counts, own original annotation vs. this pipeline's output ---
    fig, ax = plt.subplots(figsize=figsize)
    width = 0.35
    ax.bar([i - width / 2 for i in x], [r["n_original"] for r in rows], width, label="own original annotation")
    ax.bar([i + width / 2 for i in x], [r["n_final"] for r in rows], width, label="pipeline final annotation")
    ax.set_xticks(list(x))
    ax.set_xticklabels(strains, rotation=45, ha="right")
    ax.set_ylabel("Protein-coding genes")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "01_gene_counts.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 2: valid-ORF rate ---
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(strains, [r["pct_valid_orf"] for r in rows], color="#3B7DD8")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("% of genes with a valid ORF")
    ax.set_ylim(0, 100)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "02_valid_orf_rate.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 3: extra-copy genes, split by valid ORF ---
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(strains, [r["n_extra_copy_valid_orf"] for r in rows], label="valid ORF", color="#2F7A55")
    ax.bar(strains, [r["n_extra_copy_invalid_orf"] for r in rows],
           bottom=[r["n_extra_copy_valid_orf"] for r in rows], label="no valid ORF", color="#C0392B")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("Extra-copy genes")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "03_extra_copies.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 4: overlap severity ---
    fig, ax = plt.subplots(figsize=figsize)
    sev_colors = {"minor": "#B0B0B0", "substantial": "#F2B035", "contained": "#C0392B"}
    bottoms = [0] * len(strains)
    for sev in ("minor", "substantial", "contained"):
        values = [severity_all[s][sev] for s in strains]
        ax.bar(strains, values, bottom=bottoms, label=sev, color=sev_colors[sev])
        bottoms = [b + v for b, v in zip(bottoms, values)]
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("CDS-overlapping pairs")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "04_overlap_severity.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 5: source composition ---
    src_labels = ["main", "original", "reference_rescued", "good_others", "rescued_others"]
    src_colors = {"main": "#F8766D", "original": "#00BA38", "reference_rescued": "#619CFF",
                  "good_others": "#B379E0", "rescued_others": "#7A5195"}
    fig, ax = plt.subplots(figsize=figsize)
    bottoms = [0] * len(strains)
    for src in src_labels:
        values = [source_all[s][src] for s in strains]
        ax.bar(strains, values, bottom=bottoms, label=src, color=src_colors[src])
        bottoms = [b + v for b, v in zip(bottoms, values)]
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("Protein-coding genes")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "05_source_composition.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 6: overlap pairs by strand (same pairs as plot 4, split differently) ---
    fig, ax = plt.subplots(figsize=figsize)
    strand_colors = {"same_strand": "#C0392B", "opposite_strand": "#B0B0B0"}
    bottoms = [0] * len(strains)
    for key, label in (("same_strand", "same strand"), ("opposite_strand", "opposite strand")):
        values = [severity_all[s][key] for s in strains]
        ax.bar(strains, values, bottom=bottoms, label=label, color=strand_colors[key])
        bottoms = [b + v for b, v in zip(bottoms, values)]
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("CDS-overlapping pairs")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "06_overlap_strand.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 7: overlap pairs by domain support (same pairs as plot 4, split differently) ---
    fig, ax = plt.subplots(figsize=figsize)
    domain_colors = {"n_both_domain": "#2F7A55", "n_one_domain": "#F2B035", "n_no_domain": "#B0B0B0"}
    domain_labels = {"n_both_domain": "domain on both", "n_one_domain": "domain on one", "n_no_domain": "domain on neither"}
    domain_by_strain = {r["strain"]: r for r in domain_summary_rows}
    bottoms = [0] * len(strains)
    for key in ("n_both_domain", "n_one_domain", "n_no_domain"):
        values = [domain_by_strain[s][key] for s in strains]
        ax.bar(strains, values, bottom=bottoms, label=domain_labels[key], color=domain_colors[key])
        bottoms = [b + v for b, v in zip(bottoms, values)]
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("CDS-overlapping pairs")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "07_overlap_domain_support.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Plot 8: minor-severity overlaps, by position in the bigger gene's own exon structure ---
    fig, ax = plt.subplots(figsize=figsize)
    pos_colors = {"flank": "#3B7DD8", "inside": "#C0392B", "single_exon": "#B0B0B0"}
    pos_labels = {"flank": "flank (likely both real)", "inside": "inside (likely smaller gene is spurious)",
                  "single_exon": "single-exon bigger gene (n/a)"}
    bottoms = [0] * len(strains)
    for key in ("flank", "inside", "single_exon"):
        values = [minor_position_all[s][key] for s in strains]
        ax.bar(strains, values, bottom=bottoms, label=pos_labels[key], color=pos_colors[key])
        bottoms = [b + v for b, v in zip(bottoms, values)]
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_ylabel("Minor-severity CDS-overlapping pairs")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=1, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "08_minor_overlap_position.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"\nDone - 8 plots + summary.tsv + extra_copy_genes.tsv ({len(extra_copy_rows)} rows) + "
          f"cds_overlapping_pairs.tsv ({len(pair_rows)} rows) + overlap_domain_summary.tsv "
          f"written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
