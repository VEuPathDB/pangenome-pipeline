#!/usr/bin/env python3
"""
Step 3: build a sequential panlocus catalogue.

This script keeps validation evidence separate from final panlocus calls and
records native annotation structure events such as split and fused gene models.
"""

import argparse
import csv
import glob
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

from panloci_common import (
    extract_sequence,
    find_gff_file,
    is_pseudogene_feature,
    load_fasta,
    merge_hal_liftover,
    parse_gff_attributes,
    parse_samples,
    split_parent_ids,
)


STATUS_VALUES = {"PRESENT", "PARTIAL", "ABSENT"}


def extract_gene_pseudogene_bed(gff_path, output_bed, audit_tsv=None):
    """
    Extract native annotation intervals from a GFF.

    Protein-coding gene-like features are represented by the CDS span of the
    largest CDS transcript. Pseudogenes are represented by their annotated
    feature span.
    """
    gene_records = {}
    pseudogene_records = {}
    transcript_to_gene = {}
    cds_by_gene_transcript = defaultdict(lambda: defaultdict(list))

    with open(gff_path, "r", encoding="utf-8") as inp:
        for line in inp:
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


def read_bed(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path, "r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            try:
                rows.append({
                    "chrom": cols[0],
                    "start": int(cols[1]),
                    "end": int(cols[2]),
                    "name": cols[3],
                    "score": cols[4],
                    "strand": cols[5],
                    "orig_size": int(cols[6]) if len(cols) > 6 and cols[6] else int(cols[2]) - int(cols[1]),
                    "lifted_size": int(cols[7]) if len(cols) > 7 and cols[7] else int(cols[2]) - int(cols[1]),
                    "source_file": os.path.basename(path),
                    "line_no": line_no,
                })
            except ValueError:
                continue
    return rows


def annotation_bed_path(beddir, strain_id):
    return os.path.join(beddir, f"{strain_id}_cds_span_gene_pseudogene.bed")


def find_existing_annotation_bed(beddir, strain_id):
    preferred = annotation_bed_path(beddir, strain_id)
    if os.path.exists(preferred):
        return preferred
    return preferred


def write_bed(path, rows):
    with open(path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['name']}\t0\t{row['strand']}\n")


def classify_validation(
    coverage,
    identity,
    present_cov_cutoff=0.50,
    present_id_cutoff=50.0,
    partial_cov_cutoff=0.30,
    partial_id_cutoff=30.0,
):
    if coverage >= present_cov_cutoff and identity >= present_id_cutoff:
        return "PRESENT"
    if coverage >= partial_cov_cutoff and identity >= partial_id_cutoff:
        return "PARTIAL"
    return "ABSENT"


def run_minimap2_pairwise(seq_query, seq_target, locus_id, minimap_bin="minimap2", cactus_img=None):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as qf, \
            tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tf:
        qf.write(f">{locus_id}\n{seq_query}\n")
        tf.write(f">{locus_id}\n{seq_target}\n")
        q_name = qf.name
        t_name = tf.name
    try:
        if cactus_img:
            cmd = ["apptainer", "exec", cactus_img, "minimap2", "-c", t_name, q_name]
        else:
            cmd = [minimap_bin, "-c", t_name, q_name]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        best = None
        for line in res.stdout.splitlines():
            cols = line.split("\t")
            if len(cols) < 11:
                continue
            q_len = int(cols[1])
            q_start = int(cols[2])
            q_end = int(cols[3])
            matches = int(cols[9])
            block_len = int(cols[10])
            if block_len <= 0 or q_len <= 0:
                continue
            record = {
                "query_length": q_len,
                "target_length": int(cols[6]),
                "alignment_length": block_len,
                "query_coverage": round((q_end - q_start) / q_len, 4),
                "percent_identity": round((matches / block_len) * 100, 2),
                "matches": matches,
            }
            if best is None or record["matches"] > best["matches"]:
                best = record
        return best
    except subprocess.CalledProcessError as exc:
        print(f"  [Warning] minimap2 failed for {locus_id}: {exc.stderr.strip()}", file=sys.stderr)
        return None
    finally:
        for path in (q_name, t_name):
            if os.path.exists(path):
                os.remove(path)


def overlap_bp(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def build_interval_index(intervals):
    index = defaultdict(list)
    for item in intervals:
        index[(item["chrom"], item["strand"])].append(item)
    for key in index:
        index[key].sort(key=lambda item: (item["start"], item["end"]))
    return index


def find_relationships(native_genes, existing_intervals, cutoff):
    index = build_interval_index(existing_intervals)
    by_native = defaultdict(list)
    by_locus = defaultdict(list)
    for native in native_genes:
        n_len = max(1, native["end"] - native["start"])
        for hit in index.get((native["chrom"], native["strand"]), []):
            if hit["start"] >= native["end"]:
                break
            if hit["end"] <= native["start"]:
                continue
            ov = overlap_bp(native["start"], native["end"], hit["start"], hit["end"])
            if ov <= 0:
                continue
            t_len = max(1, hit["end"] - hit["start"])
            frac_native = ov / n_len
            frac_transfer = ov / t_len
            if frac_native >= cutoff or frac_transfer >= cutoff:
                rel = {
                    "native_gene_id": native["name"],
                    "native": native,
                    "panlocus_id": hit["panlocus_id"],
                    "interval": hit,
                    "overlap_bp": ov,
                    "overlap_fraction_of_native": round(frac_native, 6),
                    "overlap_fraction_of_transfer": round(frac_transfer, 6),
                }
                by_native[native["name"]].append(rel)
                by_locus[hit["panlocus_id"]].append(rel)
    return by_native, by_locus


def union_coverage_on_transfer(relationships):
    by_interval = defaultdict(list)
    lengths = {}
    for rel in relationships:
        interval = rel["interval"]
        key = (interval["chrom"], interval["start"], interval["end"], interval["strand"])
        by_interval[key].append((rel["native"]["start"], rel["native"]["end"]))
        lengths[key] = max(1, interval["end"] - interval["start"])
    best = 0.0
    for key, spans in by_interval.items():
        spans = sorted((max(start, key[1]), min(end, key[2])) for start, end in spans)
        merged = []
        for start, end in spans:
            if start >= end:
                continue
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        covered = sum(end - start for start, end in merged)
        best = max(best, covered / lengths[key])
    return round(best, 6)


def new_empty_evidence(final_call="ABSENT"):
    return {
        "final_call": final_call,
        "validation_status": "ABSENT",
        "query_coverage": "0.0",
        "percent_identity": "0.0",
        "alignment_length": "0",
        "target_length": "0",
        "evidence_type": "no_evidence",
        "confidence_class": "none",
        "review_flag": "FALSE",
        "review_reason": "NA",
        "annotation_structure": "none",
        "native_annotation_available": "FALSE",
        "native_gene_ids": "NA",
        "native_gene_count": "0",
        "matched_locus_ids": "NA",
        "best_native_overlap_fraction": "NA",
        "best_transfer_overlap_fraction": "NA",
        "combined_native_overlap_fraction": "NA",
        "coordinates": "NA",
        "source_gene_id": "NA",
        "source_strain": "NA",
    }


def coord_string(coord):
    if not coord:
        return "NA"
    return f"{coord['chrom']}:{coord['start']}-{coord['end']}({coord['strand']})"


def load_validation_table(path):
    records = {}
    if not os.path.exists(path):
        return records
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            records[row["gene_id"]] = row
    return records


def choose_best_relationship(rels):
    return sorted(
        rels,
        key=lambda rel: (
            rel["overlap_fraction_of_transfer"],
            rel["overlap_bp"],
            rel["overlap_fraction_of_native"],
        ),
        reverse=True,
    )[0]


def decide_final_from_validation_and_native(validation_status, has_native):
    if validation_status == "PRESENT":
        return "PRESENT", "transfer_validation_plus_native_annotation" if has_native else "transfer_validation", "high", "FALSE", "NA"
    if validation_status == "PARTIAL":
        if has_native:
            return "PRESENT", "transfer_partial_plus_native_annotation", "medium", "FALSE", "NA"
        return "PARTIAL", "transfer_partial_only", "medium", "FALSE", "NA"
    if has_native:
        return "PRESENT", "native_annotation_overlap_without_transfer_validation", "low", "TRUE", "native_gene_overlaps_absent_transfer"
    return "ABSENT", "no_valid_transfer", "none", "FALSE", "NA"


def validate_novel_transfer(locus_id, donor_coord, target_coord, donor_fasta, target_fasta, args):
    seq_query = extract_sequence(donor_fasta, donor_coord["chrom"], donor_coord["start"], donor_coord["end"], donor_coord["strand"])
    seq_target = extract_sequence(target_fasta, target_coord["chrom"], target_coord["start"], target_coord["end"], target_coord["strand"])
    if not seq_query or not seq_target:
        return None
    res = run_minimap2_pairwise(seq_query, seq_target, locus_id, args.minimap_bin, args.cactus_img)
    if not res:
        return None
    cov = round(target_coord["lifted_size"] / target_coord["orig_size"], 4) if target_coord.get("orig_size", 0) else 0.0
    status = classify_validation(
        cov,
        res["percent_identity"],
        args.present_cov_cutoff,
        args.present_id_cutoff,
        args.partial_cov_cutoff,
        args.partial_id_cutoff,
    )
    res["query_coverage"] = cov
    res["validation_status"] = status
    return res


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description="Step 3: build sequential panlocus catalogue")
    parser.add_argument("-s", "--samples", default="samples.tsv", help="Path to samples.tsv")
    parser.add_argument("-b", "--beddir", default="./candidate_transfers", help="Directory with candidate transfer BED files")
    parser.add_argument("-v", "--valdir", default="./validation", help="Directory with step 2 validation TSV files")
    parser.add_argument("-o", "--outdir", default="./panlocus_catalogue", help="Output directory")
    parser.add_argument("-l", "--hal", help="HAL alignment file for transferring novel panloci")
    parser.add_argument("-c", "--cactus-img", help="Cactus Apptainer/Singularity image")
    parser.add_argument("-m", "--minimap-bin", default="minimap2", help="minimap2 executable")
    parser.add_argument("-g", "--gap", type=int, default=5000, help="Maximum gap for merging liftover blocks")
    parser.add_argument("--native-overlap-cutoff", type=float, default=0.50, help="Native/transfer overlap cutoff")
    parser.add_argument("--present-cov-cutoff", type=float, default=0.50, help="Coverage cutoff for PRESENT")
    parser.add_argument("--present-id-cutoff", type=float, default=50.0, help="Identity cutoff for PRESENT")
    parser.add_argument("--partial-cov-cutoff", type=float, default=0.30, help="Coverage cutoff for PARTIAL")
    parser.add_argument("--partial-id-cutoff", type=float, default=30.0, help="Identity cutoff for PARTIAL")
    args = parser.parse_args()

    samples = parse_samples(args.samples)
    ref = samples[0]
    ref_id = ref["strain_id"]
    sample_by_id = {sample["strain_id"]: sample for sample in samples}
    order_by_id = {sample["strain_id"]: sample["analysis_order"] for sample in samples}
    strain_ids = [sample["strain_id"] for sample in samples]

    os.makedirs(args.outdir, exist_ok=True)
    occupancy_dir = os.path.join(args.outdir, "occupancy")
    reports_dir = os.path.join(args.outdir, "iteration_reports")
    qc_dir = os.path.join(args.outdir, "qc")
    os.makedirs(occupancy_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)
    os.makedirs(qc_dir, exist_ok=True)
    os.makedirs(os.path.join(args.beddir, "raw_liftovers"), exist_ok=True)

    print("======================================================================")
    print(" Step 3: Build Sequential Panlocus Catalogue")
    print("======================================================================")
    print(f"  Seed strain:                 {ref_id}")
    print(f"  Native overlap cutoff:       {args.native_overlap_cutoff}")
    print(f"  PRESENT: coverage >= {args.present_cov_cutoff}, identity >= {args.present_id_cutoff}")
    print(f"  PARTIAL: coverage >= {args.partial_cov_cutoff}, identity >= {args.partial_id_cutoff}")
    print("======================================================================")

    seed_beds = glob.glob(os.path.join(args.beddir, "*_seed.bed"))
    existing_ref_annotation_bed = find_existing_annotation_bed(args.beddir, ref_id)
    if os.path.exists(existing_ref_annotation_bed):
        seed_beds.append(existing_ref_annotation_bed)
    if not seed_beds:
        ref_gff = find_gff_file(ref["annotation_gff"], ref_id)
        if not ref_gff:
            raise RuntimeError(f"Could not find seed BED or reference GFF for {ref_id}")
        seed_bed = annotation_bed_path(args.beddir, ref_id)
        seed_audit = os.path.join(args.beddir, "qc", f"{ref_id}_cds_span_annotation_audit.tsv")
        extract_gene_pseudogene_bed(ref_gff, seed_bed, seed_audit)
    else:
        seed_bed = seed_beds[0]

    panloci = {}
    coords = defaultdict(dict)
    evidence = {}
    annotation_events = []
    review_rows = []
    novel_rows = []
    fasta_cache = {}

    for row in read_bed(seed_bed):
        locus_id = row["name"]
        panloci[locus_id] = {"panlocus_id": locus_id, "source_strain": ref_id, "source_gene_id": locus_id}
        coords[locus_id][ref_id] = row

    for locus_id in panloci:
        ev = new_empty_evidence("PRESENT")
        ev.update({
            "validation_status": "PRESENT",
            "query_coverage": "1.0",
            "percent_identity": "100.0",
            "alignment_length": str(coords[locus_id][ref_id]["end"] - coords[locus_id][ref_id]["start"]),
            "evidence_type": "seed_reference_annotation",
            "confidence_class": "high",
            "native_annotation_available": "TRUE",
            "native_gene_ids": locus_id,
            "native_gene_count": "1",
            "annotation_structure": "seed_model",
            "coordinates": coord_string(coords[locus_id][ref_id]),
            "source_gene_id": locus_id,
            "source_strain": ref_id,
        })
        evidence[(locus_id, ref_id)] = ev

    print(f"[+] Loaded {len(panloci)} seed panloci from {seed_bed}")

    # Load initial seed validation evidence and coordinates for every strain.
    for sample in samples[1:]:
        strain_id = sample["strain_id"]
        val_path = os.path.join(args.valdir, f"{ref_id}_to_{strain_id}.validation.tsv")
        val_records = load_validation_table(val_path)
        transfer_by_locus = {row["name"]: row for row in read_bed(os.path.join(args.beddir, f"{ref_id}_to_{strain_id}_merged.bed"))}
        for locus_id in panloci:
            val = val_records.get(locus_id)
            status = val["validation_status"] if val else "ABSENT"
            ev = new_empty_evidence(status)
            if val:
                ev.update({
                    "validation_status": status,
                    "query_coverage": val.get("query_coverage", "0.0"),
                    "percent_identity": val.get("percent_identity", "0.0"),
                    "alignment_length": val.get("alignment_length", "0"),
                    "target_length": val.get("target_length", "0"),
                    "evidence_type": "transfer_validation" if status != "ABSENT" else "no_valid_transfer",
                    "confidence_class": "high" if status == "PRESENT" else ("medium" if status == "PARTIAL" else "none"),
                    "source_gene_id": locus_id,
                    "source_strain": ref_id,
                })
            if status in {"PRESENT", "PARTIAL"} and locus_id in transfer_by_locus:
                coords[locus_id][strain_id] = transfer_by_locus[locus_id]
                ev["coordinates"] = coord_string(transfer_by_locus[locus_id])
            evidence[(locus_id, strain_id)] = ev

    annotated_targets = [sample for sample in samples[1:] if sample["annotation_gff"] != "NA"]
    for target in annotated_targets:
        target_id = target["strain_id"]
        print(f"\n[+] Sequential annotation integration: {target_id}")
        gff = find_gff_file(target["annotation_gff"], target_id)
        if not gff:
            print(f"  [Warning] No GFF found for {target_id}; skipping native annotation integration.", file=sys.stderr)
            continue
        native_bed = find_existing_annotation_bed(args.beddir, target_id)
        if not os.path.exists(native_bed):
            native_audit = os.path.join(args.beddir, "qc", f"{target_id}_cds_span_annotation_audit.tsv")
            extract_gene_pseudogene_bed(gff, native_bed, native_audit)
        native_genes = read_bed(native_bed)

        existing_intervals = []
        seed_transfer = os.path.join(args.beddir, f"{ref_id}_to_{target_id}_merged.bed")
        for row in read_bed(seed_transfer):
            if row["name"] in panloci:
                row["panlocus_id"] = row["name"]
                existing_intervals.append(row)
        for path in glob.glob(os.path.join(args.beddir, f"*_to_{target_id}_novel_merged.bed")):
            basename = os.path.basename(path)
            source = basename.split(f"_to_{target_id}_novel_merged.bed")[0]
            if source == target_id or order_by_id.get(source, 10**9) >= target["analysis_order"]:
                continue
            for row in read_bed(path):
                if row["name"] in panloci:
                    row["panlocus_id"] = row["name"]
                    existing_intervals.append(row)

        by_native, by_locus = find_relationships(native_genes, existing_intervals, args.native_overlap_cutoff)
        represented_native_ids = set(by_native)

        for locus_id, rels in by_locus.items():
            native_ids = sorted({rel["native_gene_id"] for rel in rels})
            matched_loci_by_native = set()
            for native_id in native_ids:
                matched_loci_by_native.update(rel["panlocus_id"] for rel in by_native[native_id])
            if len(native_ids) > 1 and len(matched_loci_by_native) > 1:
                annotation_structure = "complex_split_and_fused_model"
            elif len(native_ids) > 1:
                annotation_structure = "split_native_model"
            elif len(matched_loci_by_native) > 1:
                annotation_structure = "fused_native_model"
            else:
                annotation_structure = "single_native_model"

            validation_status = evidence.get((locus_id, target_id), new_empty_evidence())["validation_status"]
            final_call, evidence_type, confidence, review_flag, review_reason = decide_final_from_validation_and_native(validation_status, True)
            best = choose_best_relationship(rels)
            coords[locus_id][target_id] = best["native"].copy()
            coords[locus_id][target_id]["gene_id"] = best["native_gene_id"]

            current_ev = evidence.get((locus_id, target_id), new_empty_evidence())
            current_ev.update({
                "final_call": final_call,
                "evidence_type": evidence_type,
                "confidence_class": confidence,
                "review_flag": review_flag,
                "review_reason": review_reason,
                "annotation_structure": annotation_structure,
                "native_annotation_available": "TRUE",
                "native_gene_ids": ",".join(native_ids),
                "native_gene_count": str(len(native_ids)),
                "matched_locus_ids": ",".join(sorted(matched_loci_by_native)),
                "best_native_overlap_fraction": str(best["overlap_fraction_of_native"]),
                "best_transfer_overlap_fraction": str(best["overlap_fraction_of_transfer"]),
                "combined_native_overlap_fraction": str(union_coverage_on_transfer(rels)),
                "coordinates": coord_string(coords[locus_id][target_id]),
                "source_gene_id": panloci[locus_id]["source_gene_id"],
                "source_strain": panloci[locus_id]["source_strain"],
            })
            evidence[(locus_id, target_id)] = current_ev

            if annotation_structure != "single_native_model":
                annotation_events.append({
                    "strain_id": target_id,
                    "panlocus_id": locus_id,
                    "annotation_structure": annotation_structure,
                    "native_gene_ids": ",".join(native_ids),
                    "matched_locus_ids": ",".join(sorted(matched_loci_by_native)),
                    "combined_native_overlap_fraction": current_ev["combined_native_overlap_fraction"],
                    "final_call": final_call,
                    "validation_status": validation_status,
                    "review_flag": review_flag,
                })
            if review_flag == "TRUE":
                review_rows.append({
                    "strain_id": target_id,
                    "panlocus_id": locus_id,
                    "review_reason": review_reason,
                    "validation_status": validation_status,
                    "final_call": final_call,
                    "native_gene_ids": ",".join(native_ids),
                    "annotation_structure": annotation_structure,
                })

        novel_candidates = [gene for gene in native_genes if gene["name"] not in represented_native_ids]
        novel_bed = os.path.join(reports_dir, f"donor_{target_id}_novel_candidates.bed")
        write_bed(novel_bed, novel_candidates)
        print(f"  Native genes: {len(native_genes)} | represented: {len(represented_native_ids)} | novel: {len(novel_candidates)}")

        novel_id_by_native = {}
        for gene in novel_candidates:
            locus_id = gene["name"]
            if locus_id in panloci:
                locus_id = f"{target_id}_{gene['name']}"
            novel_id_by_native[gene["name"]] = locus_id
            panloci[locus_id] = {"panlocus_id": locus_id, "source_strain": target_id, "source_gene_id": gene["name"]}
            coords[locus_id][target_id] = gene
            ev = new_empty_evidence("PRESENT")
            ev.update({
                "validation_status": "PRESENT",
                "query_coverage": "1.0",
                "percent_identity": "100.0",
                "evidence_type": "native_annotation_novel_locus",
                "confidence_class": "high",
                "annotation_structure": "single_native_model",
                "native_annotation_available": "TRUE",
                "native_gene_ids": gene["name"],
                "native_gene_count": "1",
                "coordinates": coord_string(gene),
                "source_gene_id": gene["name"],
                "source_strain": target_id,
            })
            evidence[(locus_id, target_id)] = ev
            novel_rows.append({
                "source_strain": target_id,
                "panlocus_id": locus_id,
                "native_gene_id": gene["name"],
                "chrom": gene["chrom"],
                "start": gene["start"],
                "end": gene["end"],
                "strand": gene["strand"],
                "reason": "no_existing_panlocus_overlap_at_cutoff",
            })

        if novel_candidates and args.hal:
            donor_fasta = fasta_cache.get(target["genome_fasta"])
            if donor_fasta is None:
                donor_fasta = load_fasta(target["genome_fasta"])
                fasta_cache[target["genome_fasta"]] = donor_fasta
            qc_multihit = os.path.join(qc_dir, "novel_multihit_loci.tsv")
            for sample in samples:
                s_id = sample["strain_id"]
                if s_id == target_id:
                    continue
                raw_bed = os.path.join(args.beddir, "raw_liftovers", f"{target_id}_to_{s_id}_novel_raw.bed")
                merged_bed = os.path.join(args.beddir, f"{target_id}_to_{s_id}_novel_merged.bed")
                cmd = ["halLiftover", args.hal, "--bedType", "6", target_id, novel_bed, s_id, raw_bed]
                if args.cactus_img:
                    cmd = ["apptainer", "exec", args.cactus_img] + cmd
                try:
                    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
                    merge_hal_liftover(raw_bed, novel_bed, merged_bed, qc_multihit, args.gap, id_column="panlocus_id")
                except subprocess.CalledProcessError as exc:
                    print(f"  [Warning] Novel transfer failed {target_id} -> {s_id}: {exc.stderr.strip()}", file=sys.stderr)
                    continue
                target_fasta = fasta_cache.get(sample["genome_fasta"])
                if target_fasta is None:
                    target_fasta = load_fasta(sample["genome_fasta"])
                    fasta_cache[sample["genome_fasta"]] = target_fasta
                transfers_by_locus = {row["name"]: row for row in read_bed(merged_bed)}
                for novel in novel_candidates:
                    locus_id = novel_id_by_native[novel["name"]]
                    transfer = transfers_by_locus.get(novel["name"])
                    if not transfer:
                        evidence[(locus_id, s_id)] = new_empty_evidence("ABSENT")
                        continue
                    result = validate_novel_transfer(locus_id, coords[locus_id][target_id], transfer, donor_fasta, target_fasta, args)
                    if not result:
                        evidence[(locus_id, s_id)] = new_empty_evidence("ABSENT")
                        continue
                    status = result["validation_status"]
                    ev = new_empty_evidence(status)
                    ev.update({
                        "validation_status": status,
                        "query_coverage": str(result["query_coverage"]),
                        "percent_identity": str(result["percent_identity"]),
                        "alignment_length": str(result["alignment_length"]),
                        "target_length": str(result["target_length"]),
                        "evidence_type": "novel_transfer_validation" if status != "ABSENT" else "no_valid_transfer",
                        "confidence_class": "high" if status == "PRESENT" else ("medium" if status == "PARTIAL" else "none"),
                        "source_gene_id": panloci[locus_id]["source_gene_id"],
                        "source_strain": panloci[locus_id]["source_strain"],
                    })
                    if status in {"PRESENT", "PARTIAL"}:
                        coords[locus_id][s_id] = transfer
                        ev["coordinates"] = coord_string(transfer)
                    evidence[(locus_id, s_id)] = ev

    # Mark strains without annotation explicitly for any missing evidence.
    for locus_id in panloci:
        for sample in samples:
            key = (locus_id, sample["strain_id"])
            if key not in evidence:
                evidence[key] = new_empty_evidence("ABSENT")
            if sample["annotation_gff"] == "NA":
                evidence[key]["native_annotation_available"] = "FALSE"
            evidence[key]["coordinates"] = coord_string(coords.get(locus_id, {}).get(sample["strain_id"])) if coords.get(locus_id, {}).get(sample["strain_id"]) else evidence[key]["coordinates"]
            evidence[key]["source_gene_id"] = panloci[locus_id]["source_gene_id"]
            evidence[key]["source_strain"] = panloci[locus_id]["source_strain"]

    catalogue_path = os.path.join(args.outdir, "panlocus_catalogue.tsv")
    with open(catalogue_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["panlocus_id", "source_strain", "source_gene_id", "mapped_strains_count"])
        for locus_id in sorted(panloci):
            mapped = sum(1 for strain_id in strain_ids if evidence[(locus_id, strain_id)]["final_call"] in {"PRESENT", "PARTIAL"})
            writer.writerow([locus_id, panloci[locus_id]["source_strain"], panloci[locus_id]["source_gene_id"], mapped])

    long_fields = [
        "panlocus_id", "strain_id", "final_call", "validation_status", "query_coverage",
        "percent_identity", "alignment_length", "target_length", "coordinates", "evidence_type",
        "confidence_class", "review_flag", "review_reason", "annotation_structure",
        "native_annotation_available", "native_gene_ids", "native_gene_count", "matched_locus_ids",
        "best_native_overlap_fraction", "best_transfer_overlap_fraction",
        "combined_native_overlap_fraction", "source_gene_id", "source_strain",
    ]
    long_rows = []
    for locus_id in sorted(panloci):
        for strain_id in strain_ids:
            ev = evidence[(locus_id, strain_id)]
            row = {"panlocus_id": locus_id, "strain_id": strain_id}
            row.update(ev)
            long_rows.append({field: row.get(field, "NA") for field in long_fields})
    write_tsv(os.path.join(occupancy_dir, "final_occupancy_long.tsv"), long_rows, long_fields)

    with open(os.path.join(occupancy_dir, "final_occupancy_matrix.tsv"), "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["panlocus_id"] + strain_ids)
        for locus_id in sorted(panloci):
            writer.writerow([locus_id] + [evidence[(locus_id, strain_id)]["final_call"] for strain_id in strain_ids])

    write_tsv(
        os.path.join(args.outdir, "annotation_structure_events.tsv"),
        annotation_events,
        ["strain_id", "panlocus_id", "annotation_structure", "native_gene_ids", "matched_locus_ids", "combined_native_overlap_fraction", "final_call", "validation_status", "review_flag"],
    )
    write_tsv(
        os.path.join(args.outdir, "review_flags.tsv"),
        review_rows,
        ["strain_id", "panlocus_id", "review_reason", "validation_status", "final_call", "native_gene_ids", "annotation_structure"],
    )
    write_tsv(
        os.path.join(args.outdir, "novel_candidates.tsv"),
        novel_rows,
        ["source_strain", "panlocus_id", "native_gene_id", "chrom", "start", "end", "strand", "reason"],
    )

    print("\n======================================================================")
    print(" Step 3 complete")
    print(f"  Panloci: {len(panloci)}")
    print(f"  Catalogue: {catalogue_path}")
    print(f"  Long table: {os.path.join(occupancy_dir, 'final_occupancy_long.tsv')}")
    print("======================================================================")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
