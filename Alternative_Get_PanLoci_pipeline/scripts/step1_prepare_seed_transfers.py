#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
from collections import defaultdict

from panloci_common import (
    find_gff_file,
    is_pseudogene_feature,
    merge_hal_liftover,
    parse_gff_attributes,
    parse_samples,
    split_parent_ids,
)


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
            total_cds, transcript_id, span_start, span_end, exon_count = min(
                transcript_summary, key=lambda item: (-item[0], item[1])
            )
            status = "protein_coding_cds_span"
        else:
            span_start = record["start"]
            span_end = record["end"]
            transcript_id = "NA"
            total_cds = 0
            exon_count = 0
            status = "gene_without_cds_used_feature_span"
        rows.append((record["chrom"], max(0, span_start - 1), span_end, gene_id, "0", record["strand"]))
        audit_rows.append(
            [
                gene_id,
                "protein_coding_gene",
                status,
                record["chrom"],
                record["strand"],
                record["feature_type"],
                record["start"],
                record["end"],
                span_start,
                span_end,
                transcript_id,
                total_cds,
                exon_count,
                len(transcript_cds),
            ]
        )

    for gene_id, record in sorted(pseudogene_records.items()):
        rows.append((record["chrom"], max(0, record["start"] - 1), record["end"], gene_id, "0", record["strand"]))
        audit_rows.append(
            [
                gene_id,
                "pseudogene",
                "pseudogene_feature_span",
                record["chrom"],
                record["strand"],
                record["feature_type"],
                record["start"],
                record["end"],
                record["start"],
                record["end"],
                "NA",
                0,
                0,
                0,
            ]
        )

    if not rows:
        raise RuntimeError(f"No gene/pseudogene records extracted from {gff_path}")

    os.makedirs(os.path.dirname(output_bed), exist_ok=True)
    with open(output_bed, "w", encoding="utf-8") as out:
        out.writelines(
            "\t".join(map(str, row)) + "\n"
            for row in sorted(rows, key=lambda item: (item[0], item[1], item[2], item[3]))
        )

    if audit_tsv:
        os.makedirs(os.path.dirname(audit_tsv), exist_ok=True)
        with open(audit_tsv, "w", encoding="utf-8") as audit:
            audit.write(
                "gene_id\tfeature_class\tstatus\tchrom\tstrand\tfeature_type\toriginal_start\toriginal_end\tspan_start\tspan_end\tselected_transcript\ttotal_cds_length\tcds_exon_count\tn_transcripts_with_cds"
                + "\n"
            )
            audit.writelines("\t".join(map(str, row)) + "\n" for row in audit_rows)

    return len(rows)


def main():
    parser = argparse.ArgumentParser(description="Step 1: prepare seed panlocus coordinate transfers")
    parser.add_argument(
        "-s", "--samples", default="samples.tsv", help="Path to samples.tsv metadata (default: samples.tsv)"
    )
    parser.add_argument("-l", "--hal", required=True, help="Path to Cactus HAL alignment file (.hal)")
    parser.add_argument("-c", "--cactus-img", required=True, help="Path to Cactus Apptainer/Singularity image (.sif)")
    parser.add_argument(
        "-o",
        "--outdir",
        default="./candidate_transfers",
        help="Output directory for BED coordinate files (default: ./candidate_transfers)",
    )
    parser.add_argument(
        "-g",
        "--gap",
        type=int,
        default=5000,
        help="Maximum gap (bp) allowed between segments of the same gene (default: 5000 bp)",
    )
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
    samples = parse_samples(args.samples, sort_by_order=False)
    print(f"[+] Loaded {len(samples)} strain records from metadata.")
    ref_genome = None
    for s in samples:
        if s["analysis_order"] == 1:
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
        print(
            f"Error: Could not dynamically find GFF file for reference strain '{ref_id}' in directory: {ref_gff_dir}",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"  --\\> Reference GFF located: {ref_gff_path}")
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
    print(f"  --\\> Generating seed BED: {seed_bed_path}")
    try:
        audit_tsv = os.path.join(qc_dir, f"{ref_id}_cds_span_annotation_audit.tsv")
        n_seed = extract_cds_span_gene_pseudogene_bed(ref_gff_path, seed_bed_path, audit_tsv)
        print(f"  --\\> Extraction completed successfully: {n_seed} loci.")
        print(f"  --\\> Audit table: {audit_tsv}")
    except Exception as e:  # noqa: BLE001 - report any parsing failure as one line
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
        print(f"\n  --\\> Mapped Transfer: {ref_id} -> {target_strain} ...")
        raw_dir = os.path.join(args.outdir, "raw_liftovers")
        os.makedirs(raw_dir, exist_ok=True)
        raw_output_bed = os.path.join(raw_dir, f"{ref_id}_to_{target_strain}_raw.bed")
        merged_output_bed = os.path.join(args.outdir, f"{ref_id}_to_{target_strain}_merged.bed")
        cmd = [
            "apptainer",
            "exec",
            cactus_sif,
            "halLiftover",
            args.hal,
            "--bedType",
            "6",
            ref_id,
            seed_bed_path,
            target_strain,
            raw_output_bed,
        ]
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"      [SUCCESS] halLiftover: {raw_output_bed}")
            n_primary, n_dup = merge_hal_liftover(
                raw_output_bed, seed_bed_path, merged_output_bed, qc_multihit_path, args.gap
            )
            print(
                f"      [SUCCESS] Merged:    {merged_output_bed} ({n_primary} primary, {n_dup} duplicate hits logged)"
            )
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
