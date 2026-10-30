#!/usr/bin/env python3
"""Light-weight genome assembly & annotation quality assessment.

Walks an input directory containing NCBI/GenBank assemblies and computes
comparative metrics across the genomic FASTA, GFF, and protein FASTA files.
Writes two TSVs:
  1. assembly_metrics.tsv - Summary table (one row per assembly, 18 fields)
  2. protein_lengths.tsv  - Long-format table of individual protein lengths

Expected per-accession files (any GenBank/RefSeq accession):
  <acc>_*_genomic.fna      (assembled contigs/scaffolds)
  <acc>_*genomic.gff[3]    (gene / mRNA / exon / CDS annotations)
  <acc>_*protein.faa       (predicted protein sequences, optional)
  <acc>_*data_summary.tsv  (NCBI metadata summary, optional)

The pipeline layout works too: --input-dir genomes (<name>.fna) and --gff-dir gff_originals
(<name>.gff3). With no protein FASTA, the proteins are translated from the GFF's CDS.
"""

import argparse
import glob
import os
import re
import statistics
import sys
from collections import defaultdict


def parse_fasta_lengths(path):
    """Yield (header, length, gc_count, at_count) for each sequence in FASTA."""
    header = None
    length = 0
    gc = 0
    at = 0
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, length, gc, at
                header = line[1:].strip().split()[0]
                length = 0
                gc = 0
                at = 0
            else:
                seq = line.strip().upper()
                length += len(seq)
                gc += seq.count("G") + seq.count("C")
                at += seq.count("A") + seq.count("T")
        if header is not None:
            yield header, length, gc, at


def compute_n50(lengths):
    """Compute assembly N50 given a list of sequence lengths."""
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    total = sum(s)
    half = total / 2.0
    cum = 0
    for L in s:
        cum += L
        if cum >= half:
            return L
    return 0


def read_fasta(path):
    """Return {header: sequence} for a FASTA file."""
    seqs, name, parts = {}, None, []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].strip().split()[0], []
            else:
                parts.append(line.strip().upper())
    if name is not None:
        seqs[name] = "".join(parts)
    return seqs


CODON = {a + b + c: aa for (a, b, c), aa in zip(
    ((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def proteins_from_gff(gff, fna):
    """Protein lengths translated from the GFF's CDS (one per transcript), for genomes with no protein FASTA."""
    cds = defaultdict(list)
    with open(gff, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            m = re.search(r"(?:^|;)Parent=([^;]+)", fields[8])
            if m:
                phase = int(fields[7]) if fields[7] in ("0", "1", "2") else 0
                for parent in m.group(1).split(","):
                    cds[parent].append((fields[0], int(fields[3]), int(fields[4]), fields[6], phase))
    genome = read_fasta(fna)
    lengths = []
    for pieces in cds.values():
        contig, strand = pieces[0][0], pieces[0][3]
        if contig not in genome:
            continue
        pieces.sort(key=lambda p: p[1], reverse=(strand == "-"))   # 5' to 3'
        seq = "".join(genome[contig][s - 1:e] for _, s, e, _, _ in sorted(pieces, key=lambda p: p[1]))
        if strand == "-":
            seq = seq.translate(COMPLEMENT)[::-1]
        seq = seq[pieces[0][4]:]
        prot = "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
        lengths.append(len(prot[:-1] if prot.endswith("*") else prot))
    return lengths


def parse_gff(path):
    """Parse GFF3/GTF to extract feature counts and derived intron statistics.

    Introns are derived as sum(max(exons_per_mRNA - 1, 0)) across transcripts.
    Genes are the top-level features (no Parent) whose type ends in "gene" (gene, protein_coding_gene,
    ncRNA_gene, pseudogene, ...), counted as "gene" whatever the type the GFF uses.
    "CDS" counts the transcripts with a CDS (one per protein), not the CDS lines (one per coding exon).
    """
    feature_counts = defaultdict(int)
    exons_per_transcript = defaultdict(int)
    cds_transcripts = set()

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            ftype = fields[2]
            feature_counts[ftype] += 1
            if ftype != "gene" and ftype.endswith("gene") and "Parent=" not in fields[8]:
                feature_counts["gene"] += 1
            if ftype == "CDS":
                m = re.search(r"(?:^|;)\s*Parent=([^;]+)", fields[8]) or re.search(r'transcript_id "([^"]+)"', fields[8])
                if m:
                    cds_transcripts.update(m.group(1).split(","))

            if ftype == "exon":
                parent = None
                attributes = fields[8]
                for item in attributes.split(";"):
                    item = item.strip()
                    if item.startswith("Parent="):
                        parent = item.split("=", 1)[1].split(",")[0]
                        break
                    elif item.startswith("transcript_id "):
                        # GTF style support
                        parent = item.split("transcript_id ", 1)[1].strip('"').split('"')[0]
                        break
                if parent:
                    exons_per_transcript[parent] += 1

    if cds_transcripts:
        feature_counts["CDS"] = len(cds_transcripts)
    intron_count = sum(max(c - 1, 0) for c in exons_per_transcript.values())
    mean_exons = (
        statistics.mean(exons_per_transcript.values())
        if exons_per_transcript
        else (1.0 if feature_counts.get("exon", 0) > 0 else 0.0)
    )

    return feature_counts, intron_count, mean_exons


def parse_summary_tsv(path, target_acc=None):
    """Parse NCBI data_summary.tsv (single-assembly or multi-row table)."""
    if not path or not os.path.exists(path):
        return {}

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]

    if len(lines) < 2:
        return {}

    headers = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[1:]]

    # Look for matching row if multi-row
    selected_row = rows[0]
    if target_acc and len(rows) > 1:
        acc_idx = -1
        for i, h in enumerate(headers):
            if "accession" in h.lower():
                acc_idx = i
                break
        if acc_idx != -1:
            for r in rows:
                if len(r) > acc_idx and (target_acc in r[acc_idx] or r[acc_idx] in target_acc):
                    selected_row = r
                    break

    meta_raw = dict(zip(headers, selected_row))

    # Standardize common NCBI field names across datasets versions
    def get_first(candidates):
        for c in candidates:
            for k, v in meta_raw.items():
                if k.lower() == c.lower() or k.lower().replace(" ", "_") == c.lower():
                    return v
        return ""

    return {
        "organism": get_first(["Organism Scientific Name", "organism_name", "Species Name"]),
        "strain": get_first(["Organism Qualifier", "strain", "Infraspecific Name"]),
        "assembly_level": get_first(["Level", "Assembly Level", "assembly_level"]),
        "submitter_gene_count": get_first([
            "Gene Count",
            "Annotation Info Gene Count",
            "annotation_gene_count",
            "Genes",
        ]),
    }


def extract_accession(filename):
    """Extract standard GenBank/RefSeq accession (e.g. GCA_000599445.1) or basename."""
    base = os.path.basename(filename)
    m = re.match(r"^(GC[AF]_\d+\.\d+)", base)
    if m:
        return m.group(1)
    # Fallback: strip extension and common suffixes
    clean = re.sub(r"(_genomic)?\.(fna|fasta|fa|gff|gff3|faa)$", "", base)
    return clean


def discover_accessions(input_dir):
    """Discover accessions by scanning for genomic FASTA, summary TSV, or GFF."""
    accessions = set()

    patterns = [
        "*_genomic.fna",
        "*.fna",
        "*_data_summary.tsv",
        "*_genomic.gff",
        "*_genomic.gff3",
    ]
    for pat in patterns:
        for p in glob.glob(os.path.join(input_dir, pat)):
            acc = extract_accession(p)
            if acc:
                accessions.add(acc)

    return sorted(accessions)


def find_matching_file(input_dir, acc, suffixes):
    """Find the best matching file for an accession given candidate suffixes."""
    for s in suffixes:   # the exact name first, so that e.g. ME49.gff3 is not matched by ME4*
        if os.path.exists(os.path.join(input_dir, f"{acc}{s}")):
            return os.path.join(input_dir, f"{acc}{s}")
    for s in suffixes:
        matches = glob.glob(os.path.join(input_dir, f"{acc}*{s}"))
        if matches:
            return matches[0]
        # Direct exact match
        exact = os.path.join(input_dir, f"{acc}{s}")
        if os.path.exists(exact):
            return exact
    return None


def assess_one_assembly(input_dir, acc, gff_dir=None):
    """Compute complete QC metrics for a single assembly."""
    print(f"[{acc}] assessing...", file=sys.stderr)

    fna = find_matching_file(input_dir, acc, ["_genomic.fna", ".fna", ".fasta", ".fa"])
    if not fna:
        raise FileNotFoundError(f"No genomic FASTA found for accession {acc} in {input_dir}")

    gff = find_matching_file(gff_dir or input_dir, acc, ["_genomic.gff", "_genomic.gff3", ".gff3", ".gff"])
    faa = find_matching_file(input_dir, acc, ["_protein.faa", ".faa", ".faa.fasta"])

    # Summary TSV: per-accession or shared in folder
    tsv = find_matching_file(input_dir, acc, ["_data_summary.tsv", ".tsv"])
    if not tsv:
        shared_tsv = os.path.join(input_dir, "data_summary.tsv")
        if os.path.exists(shared_tsv):
            tsv = shared_tsv

    # 1. Genomic FASTA metrics
    contig_lengths = []
    gc_total = 0
    at_total = 0
    for _, L, gc, at in parse_fasta_lengths(fna):
        contig_lengths.append(L)
        gc_total += gc
        at_total += at

    genome_size = sum(contig_lengths)
    n50 = compute_n50(contig_lengths)
    denom = gc_total + at_total
    gc_pct = (100.0 * gc_total / denom) if denom > 0 else 0.0

    # 2. GFF annotations
    if gff and os.path.exists(gff):
        feature_counts, intron_count, mean_exons = parse_gff(gff)
    else:
        feature_counts, intron_count, mean_exons = defaultdict(int), 0, 0.0

    # 3. Protein lengths
    protein_lengths, protein_source = [], "none"
    if faa and os.path.exists(faa):
        protein_lengths, protein_source = [L for _, L, _, _ in parse_fasta_lengths(faa)], "protein_fasta"
    elif gff and os.path.exists(gff):
        protein_lengths, protein_source = proteins_from_gff(gff, fna), "translated_from_gff"

    protein_n = len(protein_lengths)
    protein_mean = statistics.mean(protein_lengths) if protein_lengths else 0.0
    protein_median = statistics.median(protein_lengths) if protein_lengths else 0.0

    # 4. Metadata
    meta = parse_summary_tsv(tsv, target_acc=acc) if tsv else {}

    row = {
        "accession": acc,
        "organism": meta.get("organism", ""),
        "strain": meta.get("strain", ""),
        "assembly_level": meta.get("assembly_level", ""),
        "submitter_gene_count": meta.get("submitter_gene_count", ""),
        "genome_size_bp": genome_size,
        "num_contigs": len(contig_lengths),
        "n50": n50,
        "gc_pct": round(gc_pct, 2),
        "num_genes": feature_counts.get("gene", 0),
        "num_mRNAs": feature_counts.get("mRNA", 0),
        "num_exons": feature_counts.get("exon", 0),
        "num_CDS": feature_counts.get("CDS", 0),
        "num_introns_derived": intron_count,
        "mean_exons_per_mRNA": round(mean_exons, 2),
        "num_proteins": protein_n,
        "protein_mean_len": round(protein_mean, 1),
        "protein_median_len": int(protein_median),
        "protein_source": protein_source,
    }

    return row, protein_lengths


def main():
    parser = argparse.ArgumentParser(
        description="Light-weight assessment of genome assembly & annotation quality."
    )
    parser.add_argument(
        "--input-dir",
        default=".",
        help="Input directory containing assembly FASTA, GFF, and auxiliary files (default: .)",
    )
    parser.add_argument(
        "--gff-dir",
        help="Directory with the GFFs, when they are not in --input-dir (e.g. gff_originals)",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory to write output TSV tables (default: .)",
    )
    parser.add_argument(
        "--accessions",
        nargs="*",
        help="Optional explicit list of accessions to assess (default: auto-discover all)",
    )
    args = parser.parse_args()

    if args.accessions:
        accessions = args.accessions
    else:
        accessions = discover_accessions(args.input_dir)

    if not accessions:
        print(f"Error: No assemblies found in {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Discovered {len(accessions)} assemblies in {args.input_dir}", file=sys.stderr)

    rows = []
    protein_rows = []
    for acc in accessions:
        try:
            row, lengths = assess_one_assembly(args.input_dir, acc, args.gff_dir)
            rows.append(row)
            protein_rows.extend((acc, L) for L in lengths)
        except Exception as e:
            print(f"Warning: Failed assessing {acc}: {e}", file=sys.stderr)

    if not rows:
        print("Error: No assemblies were successfully assessed.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    out_metrics = os.path.join(args.output_dir, "assembly_metrics.tsv")
    out_proteins = os.path.join(args.output_dir, "protein_lengths.tsv")

    cols = list(rows[0].keys())
    with open(out_metrics, "w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(out_proteins, "w", encoding="utf-8") as f:
        f.write("accession\tprotein_length\n")
        for acc, L in protein_rows:
            f.write(f"{acc}\t{L}\n")

    print(f"Successfully wrote {out_metrics}", file=sys.stderr)
    print(f"Successfully wrote {out_proteins}", file=sys.stderr)


if __name__ == "__main__":
    main()
