"""Helpers shared by the get_panloci step scripts (step1, step2, step3).

Not run on its own. Copy it into RUN_DIR alongside the step scripts - they import it from the
folder they are run from.
"""

import csv
import glob
import os
import sys
from collections import defaultdict

REQUIRED_SAMPLE_COLUMNS = {"strain_id", "genome_fasta", "annotation_gff", "analysis_order"}


def parse_samples(samples_tsv, sort_by_order=True):
    """samples.tsv rows as dicts: strain_id, genome_fasta, annotation_gff, analysis_order (int).

    Sorted by analysis_order unless sort_by_order is False (then file order is kept)."""
    if not os.path.exists(samples_tsv):
        raise FileNotFoundError(f"samples.tsv not found: {samples_tsv}")
    samples = []
    with open(samples_tsv, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_SAMPLE_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"samples.tsv is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            samples.append(
                {
                    "strain_id": row["strain_id"].strip(),
                    "genome_fasta": row["genome_fasta"].strip(),
                    "annotation_gff": row["annotation_gff"].strip(),
                    "analysis_order": int(row["analysis_order"]),
                }
            )
    if sort_by_order:
        samples.sort(key=lambda item: item["analysis_order"])
    return samples


def find_gff_file(path_or_dir, strain_id):
    """The GFF for a strain: path_or_dir itself if it is a file, else the *.gff/*.gff3 in that
    folder whose name contains the strain ID (or a REF-less / lower / upper-case variant of it),
    else the only GFF there. None for "NA" or when nothing matches."""
    if path_or_dir == "NA":
        return None
    if os.path.isfile(path_or_dir):
        return path_or_dir
    if not os.path.isdir(path_or_dir):
        return None
    candidates = glob.glob(os.path.join(path_or_dir, "*.gff")) + glob.glob(os.path.join(path_or_dir, "*.gff3"))
    search_terms = list(
        dict.fromkeys(
            [
                strain_id,
                strain_id.replace("REF", ""),  # e.g. REF3D7 -> 3D7
                strain_id.replace("ref", ""),
                strain_id.lower(),
                strain_id.upper(),
            ]
        )
    )
    for candidate in candidates:
        filename = os.path.basename(candidate)
        if any(term and term in filename for term in search_terms):
            return candidate
    if len(candidates) == 1:
        return candidates[0]
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
    return [item for item in parent_value.split(",") if item]


def is_pseudogene_feature(feature_type, attrs):
    feature_values = [
        feature_type,
        attrs.get("ebi_biotype", ""),
        attrs.get("biotype", ""),
        attrs.get("gene_biotype", ""),
    ]
    return any("pseudogene" in value.lower() for value in feature_values)


def load_fasta(fasta_path, missing_ok=False):
    """{first word of header: upper-case sequence}. A missing file raises FileNotFoundError, or
    with missing_ok prints an error and returns None (for callers that skip that genome)."""
    if not os.path.exists(fasta_path):
        if missing_ok:
            print(f"Error: FASTA file not found at '{fasta_path}'", file=sys.stderr)
            return None
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    fasta = {}
    name = None
    seq = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    fasta[name] = "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        fasta[name] = "".join(seq).upper()
    return fasta


def extract_sequence(fasta, chrom, start, end, strand="+"):
    """fasta[chrom][start:end] (0-based, half-open, clipped to the sequence), reverse-complemented
    for strand "-". None if the chromosome is absent or the interval is empty."""
    if not fasta or chrom not in fasta:
        return None
    chrom_len = len(fasta[chrom])
    start = max(0, min(start, chrom_len))
    end = max(0, min(end, chrom_len))
    if start >= end:
        return None
    seq = fasta[chrom][start:end]
    if strand == "-":
        comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        seq = seq.translate(comp)[::-1]
    return seq


def merge_overlapping_intervals(intervals):
    """Merge overlapping (start, end) intervals; touching half-open ones ([0,10), [10,20)) merge too."""
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


def parse_seed_lengths(seed_bed_path):
    """{name: total bp} from a BED (col 4 = name), summing a name's blocks. Missing file -> {};
    rows with fewer than 4 columns or non-numeric coordinates are skipped."""
    lengths = {}
    if not os.path.exists(seed_bed_path):
        return lengths
    with open(seed_bed_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 4:
                continue
            try:
                start, end = int(cols[1]), int(cols[2])
            except ValueError:
                continue
            lengths[cols[3]] = lengths.get(cols[3], 0) + (end - start)
    return lengths


def cluster_mapped_blocks(rows, max_gap=5000):
    """Group BED rows by (chrom, strand) and cluster their blocks: consecutive blocks at most
    max_gap bp apart (next start - previous end) share a cluster. -> [((chrom, strand), [blocks])]"""
    by_group = defaultdict(list)
    for row in rows:
        by_group[(row[0], row[5] if len(row) > 5 else "+")].append((int(row[1]), int(row[2])))
    clusters = []
    for (chrom, strand), blocks in by_group.items():
        blocks.sort()
        current = [blocks[0]]
        for block in blocks[1:]:
            if block[0] - current[-1][1] <= max_gap:
                current.append(block)
            else:
                clusters.append(((chrom, strand), current))
                current = [block]
        clusters.append(((chrom, strand), current))
    return clusters


def merge_hal_liftover(
    raw_bed_path, seed_bed_path, output_bed_path, qc_multihit_path, max_gap=5000, id_column="gene_id"
):
    """Resolve a raw halLiftover BED into one interval per locus:
      1. group the raw blocks by name (col 4);
      2. cluster blocks within max_gap bp on the same chromosome and strand;
      3. keep the cluster with the most aligned bases as the locus (output_bed_path, sorted);
      4. append every other cluster to qc_multihit_path (header's ID column named id_column).
    Returns (n_primary, n_duplicate_hits)."""
    orig_lengths = parse_seed_lengths(seed_bed_path)
    groups = defaultdict(list)
    with open(raw_bed_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                cols = line.split()
            if len(cols) < 4:
                continue
            groups[cols[3]].append(cols)

    merged = []
    duplicates = []
    for name, rows in groups.items():
        orig_size = orig_lengths.get(name, 0)
        clusters = []
        for (chrom, strand), blocks in cluster_mapped_blocks(rows, max_gap):
            mapped_len = sum(end - start for start, end in merge_overlapping_intervals(blocks))
            clusters.append(
                {
                    "chrom": chrom,
                    "strand": strand,
                    "start": min(start for start, _ in blocks),
                    "end": max(end for _, end in blocks),
                    "mapped_len": mapped_len,
                }
            )
        clusters.sort(key=lambda item: item["mapped_len"], reverse=True)
        primary = clusters[0]
        pct = round((primary["mapped_len"] / orig_size) * 100, 2) if orig_size > 0 else 0.0
        merged.append(
            (
                primary["chrom"],
                primary["start"],
                primary["end"],
                name,
                "0",
                primary["strand"],
                orig_size,
                primary["mapped_len"],
                pct,
            )
        )
        for rank, dup in enumerate(clusters[1:], start=2):
            pct_dup = round((dup["mapped_len"] / orig_size) * 100, 2) if orig_size > 0 else 0.0
            duplicates.append(
                (
                    name,
                    dup["chrom"],
                    dup["start"],
                    dup["end"],
                    dup["strand"],
                    orig_size,
                    dup["mapped_len"],
                    pct_dup,
                    rank,
                )
            )

    merged.sort(key=lambda item: (item[0], item[1]))
    with open(output_bed_path, "w", encoding="utf-8") as out:
        for row in merged:
            out.write("\t".join(map(str, row)) + "\n")
    if duplicates:
        exists = os.path.exists(qc_multihit_path)
        with open(qc_multihit_path, "a", encoding="utf-8") as out:
            if not exists:
                out.write(
                    f"{id_column}\tdup_chrom\tdup_start\tdup_end\tdup_strand\torig_size\tdup_mapped_len\tdup_percentage\thit_rank\n"
                )
            for row in duplicates:
                out.write("\t".join(map(str, row)) + "\n")
    return len(merged), len(duplicates)
