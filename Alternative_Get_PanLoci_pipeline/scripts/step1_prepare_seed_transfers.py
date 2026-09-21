#!/usr/bin/env python3
import os
import sys
import subprocess
import csv
import argparse
import glob
from collections import defaultdict

def parse_samples(samples_tsv):
    """Parses samples.tsv to extract strain information."""
    samples = []
    if not os.path.exists(samples_tsv):
        print(f"Error: Metadata file not found at '{samples_tsv}'", file=sys.stderr)
        sys.exit(1)
    with open(samples_tsv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            samples.append({
                "strain_id": row["strain_id"].strip(),
                "genome_fasta": row["genome_fasta"].strip(),
                "annotation_gff": row["annotation_gff"].strip(),
                "analysis_order": row["analysis_order"].strip()
            })
    return samples

def find_gff_file(directory, strain_id):
    """
    Dynamically finds a GFF file for the given strain inside a directory.
    Searches for common variations of the strain_id in the filename.
    """
    if not os.path.isdir(directory):
        if os.path.isfile(directory):
            return directory
        return None
    gff_files = glob.glob(os.path.join(directory, "*.gff")) + glob.glob(os.path.join(directory, "*.gff3"))
    search_terms = [
        strain_id,
        strain_id.replace('REF', ''),       # e.g., REF3D7 -> 3D7
        strain_id.replace('ref', ''),
        strain_id.lower(),
        strain_id.upper()
    ]
    search_terms = list(dict.fromkeys(search_terms))
    for file_path in gff_files:
        filename = os.path.basename(file_path)
        for term in search_terms:
            if term in filename:
                return file_path
    if len(gff_files) == 1:
        return gff_files[0]
    return None

def parse_gff_attributes(attr_text):
    """Parse a GFF3 attributes field into a dictionary."""
    attrs = {}
    for part in attr_text.strip().split(";"):
        if not part or "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key] = value
    return attrs

def split_parent_ids(parent_value):
    if not parent_value:
        return []
    return [item.strip() for item in parent_value.split(",") if item.strip()]

def is_pseudogene_feature(feature_type, attrs):
    feature_values = [
        feature_type,
        attrs.get("ebi_biotype", ""),
        attrs.get("biotype", ""),
        attrs.get("gene_biotype", ""),
    ]
    return any("pseudogene" in value.lower() for value in feature_values)

def extract_cds_span_gene_pseudogene_bed(gff_path, output_bed, audit_tsv=None):
    """
    Extract panlocus seed intervals from a GFF.

    Protein-coding gene-like features are represented by the CDS span of the
    largest CDS transcript. Pseudogenes are represented by the annotated
    pseudogene feature span.
    """
    gene_records = {}
    pseudogene_records = {}
    transcript_to_gene = {}
    cds_by_gene_transcript = defaultdict(lambda: defaultdict(list))

    with open(gff_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature_type = cols[2]
            attrs = parse_gff_attributes(cols[8])
            feature_id = attrs.get("ID") or attrs.get("Name")
            if feature_type in {"gene", "protein_coding_gene", "pseudogene"} and feature_id:
                record = {
                    "gene_id": feature_id,
                    "chrom": cols[0],
                    "start": int(cols[3]),
                    "end": int(cols[4]),
                    "strand": cols[6],
                    "feature_type": feature_type,
                }
                if is_pseudogene_feature(feature_type, attrs):
                    pseudogene_records[feature_id] = record
                else:
                    gene_records[feature_id] = record
            if feature_type in {"mRNA", "transcript"} and feature_id:
                parents = split_parent_ids(attrs.get("Parent"))
                if parents:
                    transcript_to_gene[feature_id] = parents[0]
            if feature_type == "CDS":
                parents = split_parent_ids(attrs.get("Parent"))
                gene_id = attrs.get("gene_id")
                transcript_id = parents[0] if parents else attrs.get("ID", "unknown_transcript")
                if not gene_id:
                    if transcript_id in transcript_to_gene:
                        gene_id = transcript_to_gene[transcript_id]
                    elif transcript_id in gene_records:
                        gene_id = transcript_id
                if gene_id:
                    cds_by_gene_transcript[gene_id][transcript_id].append((int(cols[3]), int(cols[4])))

    rows = []
    audit_rows = []

    for gene_id, record in sorted(gene_records.items()):
        transcript_cds = cds_by_gene_transcript.get(gene_id, {})
        if transcript_cds:
            transcript_summary = []
            for transcript_id, intervals in transcript_cds.items():
                total_cds = sum(end - start + 1 for start, end in intervals)
                span_start = min(start for start, end in intervals)
                span_end = max(end for start, end in intervals)
                transcript_summary.append((total_cds, transcript_id, span_start, span_end, len(intervals)))
            total_cds, transcript_id, span_start, span_end, exon_count = sorted(
                transcript_summary, key=lambda item: (-item[0], item[1])
            )[0]
            status = "protein_coding_cds_span"
        else:
            span_start = record["start"]
            span_end = record["end"]
            transcript_id = "NA"
            total_cds = 0
            exon_count = 0
            status = "gene_without_cds_used_feature_span"
        rows.append((record["chrom"], max(0, span_start - 1), span_end, gene_id, "0", record["strand"]))
        audit_rows.append([
            gene_id, "protein_coding_gene", status, record["chrom"], record["strand"],
            record["feature_type"], record["start"], record["end"], span_start, span_end,
            transcript_id, total_cds, exon_count, len(transcript_cds)
        ])

    for gene_id, record in sorted(pseudogene_records.items()):
        rows.append((record["chrom"], max(0, record["start"] - 1), record["end"], gene_id, "0", record["strand"]))
        audit_rows.append([
            gene_id, "pseudogene", "pseudogene_feature_span", record["chrom"], record["strand"],
            record["feature_type"], record["start"], record["end"], record["start"], record["end"],
            "NA", 0, 0, 0
        ])

    if not rows:
        raise RuntimeError(f"No gene/pseudogene records extracted from {gff_path}")

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)
    with open(output_bed, "w", encoding="utf-8") as out:
        for row in sorted(rows, key=lambda item: (item[0], item[1], item[2], item[3])):
            out.write("\t".join(map(str, row)) + "\n")

    if audit_tsv:
        os.makedirs(os.path.dirname(audit_tsv), exist_ok=True)
        with open(audit_tsv, "w", encoding="utf-8") as audit:
            audit.write("\t".join([
                "gene_id", "feature_class", "status", "chrom", "strand", "feature_type",
                "original_start", "original_end", "span_start", "span_end",
                "selected_transcript", "total_cds_length", "cds_exon_count",
                "n_transcripts_with_cds"
            ]) + "\n")
            for row in audit_rows:
                audit.write("\t".join(map(str, row)) + "\n")

    return len(rows)

def parse_seed_lengths(seed_bed_path):
    """Reads the seed BED file to build a dictionary of original gene sizes (end - start)."""
    lengths = {}
    if not os.path.exists(seed_bed_path):
        return lengths
    with open(seed_bed_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) >= 4:
                name = cols[3]
                start = int(cols[1])
                end = int(cols[2])
                lengths[name] = lengths.get(name, 0) + (end - start)
    return lengths

def cluster_mapped_blocks(rows, max_gap=5000):
    """
    Groups raw BED blocks by (chrom, strand) and clusters them.
    Two blocks are in the same cluster if they are on the same chromosome/strand
    and the distance between block1 end and block2 start is <= max_gap (default: 1000 bp).
    """
    by_group = defaultdict(list)
    for r in rows:
        chrom = r[0]
        start = int(r[1])
        end = int(r[2])
        strand = r[5] if len(r) > 5 else "+"
        by_group[(chrom, strand)].append((start, end))
    clusters = []
    for (chrom, strand), blocks in by_group.items():
        blocks.sort()
        current_cluster = [blocks[0]]
        for next_block in blocks[1:]:
            if next_block[0] - current_cluster[-1][1] <= max_gap:
                current_cluster.append(next_block)
            else:
                clusters.append(((chrom, strand), current_cluster))
                current_cluster = [next_block]
        clusters.append(((chrom, strand), current_cluster))
    return clusters

def merge_overlapping_intervals(intervals):
    """Merges overlapping or adjacent intervals and returns a list of merged (start, end) tuples."""
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda x: x[0])
    merged = [sorted_intervals[0]]
    for current in sorted_intervals[1:]:
        prev_start, prev_end = merged[-1]
        curr_start, curr_end = current
        if curr_start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, curr_end))
        else:
            merged.append(current)
    return merged

def merge_hal_liftover(raw_bed_path, seed_bed_path, output_bed_path, qc_multihit_path, max_gap=5000):
    """
    Resolves raw coordinate transfers using distance-based clustering:
      1. Groups raw blocks by gene name (ID).
      2. Clusters blocks within 5000bp on the same chromosome & strand.
      3. Selects the PRIMARY cluster (the one with the most total aligned bases) for occupancy.
      4. Saves secondary/distal clusters to a duplicate/multi-hit QC table.
    """
    # 1. Parse original lengths
    orig_lengths = parse_seed_lengths(seed_bed_path)
    # 2. Group raw liftover blocks by name
    groups = defaultdict(list)
    with open(raw_bed_path, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                cols = line.split()
            if len(cols) < 4:
                continue
            name = cols[3]
            groups[name].append(cols)
    # 3. Perform clustering and select primary hit
    merged_genes = []
    qc_duplicates = []
    for name, rows in groups.items():
        orig_size = orig_lengths.get(name, 0)
        clusters = cluster_mapped_blocks(rows, max_gap)
        processed_clusters = []
        for (chrom, strand), blocks in clusters:
            c_start = min(b[0] for b in blocks)
            c_end = max(b[1] for b in blocks)
            merged_blocks = merge_overlapping_intervals(blocks)
            c_mapped_len = sum(end - start for start, end in merged_blocks)
            processed_clusters.append({
                "chrom": chrom,
                "strand": strand,
                "start": c_start,
                "end": c_end,
                "mapped_len": c_mapped_len
            })
        processed_clusters.sort(key=lambda x: x["mapped_len"], reverse=True)
        primary = processed_clusters[0]
        percentage = round((primary["mapped_len"] / orig_size) * 100, 2) if orig_size > 0 else 0.0
        merged_genes.append((
            primary["chrom"],
            primary["start"],
            primary["end"],
            name,
            "0",
            primary["strand"],
            orig_size,
            primary["mapped_len"],
            percentage
        ))
        if len(processed_clusters) > 1:
            for rank, dup in enumerate(processed_clusters[1:], start=2):
                dup_percentage = round((dup["mapped_len"] / orig_size) * 100, 2) if orig_size > 0 else 0.0
                qc_duplicates.append((
                    name,
                    dup["chrom"],
                    dup["start"],
                    dup["end"],
                    dup["strand"],
                    orig_size,
                    dup["mapped_len"],
                    dup_percentage,
                    rank
                ))
    merged_genes.sort(key=lambda x: (x[0], x[1]))
    with open(output_bed_path, "w") as f:
        for g in merged_genes:
            f.write(f"{g[0]}\t{g[1]}\t{g[2]}\t{g[3]}\t{g[4]}\t{g[5]}\t{g[6]}\t{g[7]}\t{g[8]}\n")
    if qc_duplicates:
        file_exists = os.path.exists(qc_multihit_path)
        with open(qc_multihit_path, "a") as f:
            if not file_exists:
                f.write("gene_id\tdup_chrom\tdup_start\tdup_end\tdup_strand\torig_size\tdup_mapped_len\tdup_percentage\thit_rank\n")
            for d in qc_duplicates:
                f.write(f"{d[0]}\t{d[1]}\t{d[2]}\t{d[3]}\t{d[4]}\t{d[5]}\t{d[6]}\t{d[7]}\t{d[8]}\n")
    print(f"      [SUCCESS] Merged:    {output_bed_path} ({len(merged_genes)} primary, {len(qc_duplicates)} duplicate hits logged)")

def main():
    parser = argparse.ArgumentParser(description="Step 1: prepare seed panlocus coordinate transfers")
    parser.add_argument("-s", "--samples", default="samples.tsv", help="Path to samples.tsv metadata (default: samples.tsv)")
    parser.add_argument("-l", "--hal", required=True, help="Path to Cactus HAL alignment file (.hal)")
    parser.add_argument("-c", "--cactus-img", required=True, help="Path to Cactus Apptainer/Singularity image (.sif)")
    parser.add_argument("-o", "--outdir", default="./candidate_transfers", help="Output directory for BED coordinate files (default: ./candidate_transfers)")
    parser.add_argument("-g", "--gap", type=int, default=5000, help="Maximum gap (bp) allowed between segments of the same gene (default: 5000 bp)")
    args = parser.parse_args()
    print("======================================================================")
    print(" Step 1: Prepare Seed Panlocus Coordinate Transfers")
    print("======================================================================")
    print(f"  HAL Alignment:      {args.hal}")
    print(f"  Cactus SIF Image:   {args.cactus_img}")
    print(f"  Metadata Registry:  {args.samples}")
    print(f"  Output Directory:   {args.outdir}")
    print(f"  Maximum Merge Gap:  {args.gap} bp")
    print("======================================================================\n")
    samples = parse_samples(args.samples)
    print(f"[+] Loaded {len(samples)} strain records from metadata.")
    ref_genome = None
    for s in samples:
        if s["analysis_order"] == "1":
            ref_genome = s
            break
    if not ref_genome:
        print("Error: Could not find reference genome with analysis_order = 1 in samples.tsv", file=sys.stderr)
        sys.exit(1)
    ref_id = ref_genome["strain_id"]
    print(f"[+] Seed Genome identified: {ref_id}")
    ref_gff_dir = ref_genome["annotation_gff"]
    ref_gff_path = find_gff_file(ref_gff_dir, ref_id)
    if not ref_gff_path:
        print(f"Error: Could not dynamically find GFF file for reference strain '{ref_id}' in directory: {ref_gff_dir}", file=sys.stderr)
        sys.exit(1)
    print(f"  --\> Reference GFF located: {ref_gff_path}")
    gff_basename = os.path.basename(ref_gff_path)
    prefix = os.path.splitext(gff_basename)[0]
    seed_bed_path = os.path.join(args.outdir, f"{prefix}_seed.bed")
    os.makedirs(args.outdir, exist_ok=True)
    qc_dir = os.path.join(args.outdir, "qc")
    os.makedirs(qc_dir, exist_ok=True)
    qc_multihit_path = os.path.join(qc_dir, "multihit_duplicated_loci.tsv")
    if os.path.exists(qc_multihit_path):
        os.remove(qc_multihit_path)
    print("\n[+] Extracting seed panlocus coordinates to BED...")
    print(f"  --\> Generating seed BED: {seed_bed_path}")
    try:
        audit_tsv = os.path.join(qc_dir, f"{ref_id}_cds_span_annotation_audit.tsv")
        n_seed = extract_cds_span_gene_pseudogene_bed(ref_gff_path, seed_bed_path, audit_tsv)
        print(f"  --\> Extraction completed successfully: {n_seed} loci.")
        print(f"  --\> Audit table: {audit_tsv}")
    except Exception as e:
        print(f"Error extracting seed coordinates from GFF: {e}", file=sys.stderr)
        sys.exit(1)
    print("\n[+] Beginning coordinate transfer via halLiftover & distance-based merging...")
    success_count = 0
    fail_count = 0
    cactus_sif = args.cactus_img
    for s in samples:
        target_strain = s["strain_id"]
        if target_strain == ref_id:
            continue
        print(f"\n  --\> Mapped Transfer: {ref_id} -> {target_strain} ...")
        raw_dir = os.path.join(args.outdir, "raw_liftovers")
        os.makedirs(raw_dir, exist_ok=True)
        raw_output_bed = os.path.join(raw_dir, f"{ref_id}_to_{target_strain}_raw.bed")
        merged_output_bed = os.path.join(args.outdir, f"{ref_id}_to_{target_strain}_merged.bed")
        cmd = [
            "apptainer", "exec", cactus_sif, "halLiftover",
            args.hal, "--bedType", "6",
            ref_id, seed_bed_path,
            target_strain, raw_output_bed
        ]
        try:
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            print(f"      [SUCCESS] halLiftover: {raw_output_bed}")
            merge_hal_liftover(raw_output_bed, seed_bed_path, merged_output_bed, qc_multihit_path, args.gap)
            success_count += 1
        except subprocess.CalledProcessError as e:
            print(f"      [FAILURE] Failed to transfer. Error: {e.stderr.strip()}", file=sys.stderr)
            fail_count += 1
    print("\n======================================================================")
    print(" Seeding and Coordinate Transfer Summary")
    print("======================================================================")
    print(f"  Reference Genome: {ref_id}")
    print(f"  Seed BED created: {os.path.exists(seed_bed_path)}")
    print(f"  Successful Transfers: {success_count} / {len(samples) - 1}")
    print(f"  Failed Transfers:     {fail_count} / {len(samples) - 1}")
    print("======================================================================")

if __name__ == "__main__":
    main()
