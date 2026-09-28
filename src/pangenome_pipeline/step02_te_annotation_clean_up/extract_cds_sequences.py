#!/usr/bin/env python3
"""
Extract the spliced CDS nucleotide sequence of flagged genes from their
genome assembly FASTA, so they can be inspected directly (BLASTed, aligned
against the curated TE library, translated, etc.) rather than just read off
as coordinates.

For each flagged protein_coding or pseudogene gene, every transcript's CDS
exons are concatenated in genomic order and reverse-complemented if the gene
is on the minus strand - the standard "give me the coding sequence"
splicing, no translation. A gene with more than one transcript gets one
FASTA record per transcript; a pseudogene with no CDS feature at all (common
- see te_gene_overlap.py) contributes none, silently.

Usage:
    python3 -m pangenome_pipeline.step02_te_annotation_clean_up.extract_cds_sequences \
        --gene-table te_annotation_integration/gene_te_overlap.tsv \
        --gff-indir gff_originals \
        --genomes-dir genomes \
        --outdir te_annotation_integration/flagged_gff \
        [--min-pct 20] [--max-pct 80] [--extra-genes confirmed.tsv]

Defaults to the 20-80% "inspect" band on protein_coding and pseudogene
genes, same as extract_flagged_genes.py - --min-pct/--max-pct override that,
and --extra-genes (2-column TSV, genome_label<TAB>gene_id, no header) adds
specific genes on top, same convention as reclassify_te_genes.py.

Expects --genomes-dir to contain one assembly FASTA per genome, named
<genome_label>_<anything>.fna (e.g. GCA_000002495.2_MG8_genomic.fna for
genome_label "GCA_000002495.2") - the prefix before the first "_" block that
matches genome_label is what's looked up; contig IDs in the FASTA headers
must match the GFF's seqid column (true for genomes downloaded from NCBI
alongside their GFF3, as in this project).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pangenome_pipeline.step02_te_annotation_clean_up.reclassify_te_genes import load_extra_genes
from pangenome_pipeline.step02_te_annotation_clean_up.te_gene_overlap import index_gff_dir, parse_gff, strip_prefix

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def find_genome_fasta(genomes_dir: Path, genome_label: str) -> Path | None:
    patterns = [
        f"{genome_label}.fna",
        f"{genome_label}.fasta",
        f"{genome_label}.fa",
        f"{genome_label}_*.fna",
        f"{genome_label}_*.fasta",
        f"{genome_label}_*.fa",
        f"{genome_label}*.fna",
        f"{genome_label}*.fasta",
        f"{genome_label}*.fa",
    ]
    for pat in patterns:
        matches = sorted(genomes_dir.glob(pat))
        if matches:
            return matches[0]
    return None


def load_contigs(fasta_path: Path, wanted: set) -> dict:
    """Streams the FASTA, keeping only sequences for contigs in `wanted` -
    genome assemblies here are ~40 MB but a handful of genes only ever sit
    on a handful of contigs, so there's no need to hold the whole thing."""
    seqs = {}
    current_id = None
    current_chunks: list[str] = []
    with fasta_path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if current_id is not None and current_id in wanted:
                    seqs[current_id] = "".join(current_chunks)
                current_id = line[1:].split()[0]
                current_chunks = []
            elif current_id in wanted:
                current_chunks.append(line.strip())
        if current_id is not None and current_id in wanted:
            seqs[current_id] = "".join(current_chunks)
    return seqs


def splice_cds(contig_seq: str, cds_intervals: list, strand: str) -> str:
    """cds_intervals: list of (start, end), 1-based inclusive, any order.
    Concatenates in ascending genomic order (5'->3' for a + strand gene),
    then reverse-complements the whole thing for a - strand gene - reverse-
    complementing an ascending-order concatenation reorders and complements
    every piece correctly in one step."""
    pieces = [contig_seq[s - 1 : e] for s, e in sorted(cds_intervals)]
    seq = "".join(pieces)
    return revcomp(seq) if strand == "-" else seq


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-table", required=True, help="gene_te_overlap.tsv produced by te_gene_overlap.py")
    ap.add_argument("--gff-indir", required=True, help="Directory with the source GFF3 files (genome_label.gff/.gff3)")
    ap.add_argument(
        "--genomes-dir", required=True, help="Directory with one assembly FASTA per genome (genome_label_*.fna)"
    )
    ap.add_argument("--outdir", required=True, help="Directory to write the combined CDS FASTA to")
    ap.add_argument("--min-pct", type=float, default=20.0, help="Minimum TE_overlap_pct, inclusive (default: 20)")
    ap.add_argument("--max-pct", type=float, default=80.0, help="Maximum TE_overlap_pct, exclusive (default: 80)")
    ap.add_argument(
        "--min-div",
        type=float,
        default=None,
        help="Optional minimum TE_divergence_pct, inclusive - on top of "
        "--min-pct/--max-pct, to isolate e.g. high-overlap genes whose "
        "dominant TE fragment is divergent enough to warrant a manual look "
        "before treating the overlap call as settled",
    )
    ap.add_argument("--max-div", type=float, default=None, help="Optional maximum TE_divergence_pct, exclusive")
    ap.add_argument(
        "--extra-genes",
        default=None,
        help="Optional 2-column TSV (genome_label<TAB>gene_id, no header) of "
        "genes to add on top of the --min-pct/--max-pct band",
    )
    ap.add_argument(
        "--out-name",
        default="flagged_genes_cds.fna",
        help="Output FASTA filename, written inside --outdir (default: flagged_genes_cds.fna)",
    )
    args = ap.parse_args()

    gff_indir = Path(args.gff_indir)
    genomes_dir = Path(args.genomes_dir)
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
    sub = df[mask]
    target_pairs = set(zip(sub["genome_label"], sub["gene_id"]))

    if args.extra_genes:
        extra_pairs = load_extra_genes(args.extra_genes)
        valid_pairs = set(zip(df["genome_label"], df["gene_id"]))
        for pair in extra_pairs:
            if pair not in valid_pairs:
                print(f"WARNING: --extra-genes entry {pair} not found in {args.gene_table}, skipping")
                continue
            target_pairs.add(pair)

    if not target_pairs:
        raise SystemExit(
            f"No protein_coding/pseudogene genes with {args.min_pct}% <= "
            f"TE_overlap_pct < {args.max_pct}% found in {args.gene_table} "
            f"(plus none/no valid entries in --extra-genes)"
        )

    meta = df.set_index(["genome_label", "gene_id"])
    by_genome = {}
    for genome_label, gene_id in target_pairs:
        by_genome.setdefault(genome_label, set()).add(gene_id)

    gff_by_stem = index_gff_dir(gff_indir)
    records = []
    n_genes_written = 0
    written_pairs = set()  # (genome_label, gene_id) that got >=1 sequence written
    for genome_label, gene_ids in sorted(by_genome.items()):
        gff_path = gff_by_stem.get(genome_label)
        if gff_path is None:
            print(
                f"WARNING: no GFF found for genome_label {genome_label!r} in {gff_indir} "
                f"- skipping {len(gene_ids)} genes"
            )
            continue
        fasta_path = find_genome_fasta(genomes_dir, genome_label)
        if fasta_path is None:
            print(
                f"WARNING: no genome FASTA found for genome_label {genome_label!r} in "
                f"{genomes_dir} - skipping {len(gene_ids)} genes"
            )
            continue

        gene_meta, rna_to_gene, _, cds_by_rna, _ = parse_gff(gff_path)
        rnas_by_gene = {}
        for rna_id, parent_gene in rna_to_gene.items():
            rnas_by_gene.setdefault(strip_prefix(parent_gene, "gene-"), []).append(rna_id)

        # gene_te_overlap.tsv's gene_id is the GFF's own gene ID with a "gene-"
        # prefix stripped IF it had one (te_gene_overlap.py's strip_prefix call) -
        # NCBI GFFs prefix gene IDs with "gene-", VEuPathDB/FungiDB GFFs (e.g.
        # AF293, A1160, A1163 here) don't prefix them at all. Reconstructing the
        # raw ID by re-adding a hardcoded "gene-" only round-trips for the NCBI
        # case; build the reverse map from gene_meta itself instead, so it works
        # for either convention (or any other prefix strip_prefix is ever taught).
        raw_id_by_stripped = {strip_prefix(raw, "gene-"): raw for raw in gene_meta}

        wanted_contigs = {gene_meta[raw_id_by_stripped[gid]]["contig"] for gid in gene_ids if gid in raw_id_by_stripped}
        contig_seqs = load_contigs(fasta_path, wanted_contigs)

        for gene_id in sorted(gene_ids):
            raw_gene_id = raw_id_by_stripped.get(gene_id)
            if raw_gene_id is None:
                print(f"WARNING: gene {gene_id!r} not found in {gff_path.name} - skipping")
                continue
            g = gene_meta[raw_gene_id]
            contig_seq = contig_seqs.get(g["contig"])
            if contig_seq is None:
                print(f"WARNING: contig {g['contig']!r} not found in {fasta_path.name} for gene {gene_id!r} - skipping")
                continue
            row = meta.loc[(genome_label, gene_id)]
            transcripts = sorted(rnas_by_gene.get(gene_id, []))
            n_records_before = len(records)
            for rna_id in transcripts:
                cds_intervals = cds_by_rna.get(rna_id)
                if not cds_intervals:
                    continue
                seq = splice_cds(contig_seq, cds_intervals, g["strand"])
                transcript_label = strip_prefix(strip_prefix(rna_id, "rna-"), "mrna_")
                div = row.get("TE_divergence_pct")
                div_str = f"{div:.1f}" if pd.notna(div) else "NA"
                header = (
                    f">{genome_label}|{gene_id}|{transcript_label} "
                    f"contig={g['contig']}:{g['start']}-{g['end']}({g['strand']}) "
                    f"TE_overlap_pct={row['TE_overlap_pct']:.2f} "
                    f"TE_divergence_pct={div_str} "
                    f"dominant_TE={row['dominant_TE_class']}:{row['dominant_TE_name']} "
                    f"len={len(seq)}bp"
                )
                records.append((header, seq))
            if len(records) == n_records_before:
                # gene was found, but none of its transcripts had a CDS to splice
                # (expected for some pseudogenes - see te_gene_overlap.py) -
                # warn so a genuinely missing sequence is never just silent.
                print(
                    f"WARNING: gene {gene_id!r} ({genome_label}) has no CDS in "
                    f"{gff_path.name} across {len(transcripts)} transcript(s) - "
                    f"no sequence written"
                )
            else:
                n_genes_written += 1
                written_pairs.add((genome_label, gene_id))

    missing_pairs = sorted(target_pairs - written_pairs)
    if missing_pairs:
        print(
            f"\nWARNING: {len(missing_pairs)} of {len(target_pairs)} requested gene(s) "
            f"have NO sequence in the output FASTA (see the WARNINGs above for why - "
            f"missing GFF/FASTA, gene not found, contig not found, or no CDS):"
        )
        for genome_label, gene_id in missing_pairs:
            print(f"  {genome_label}\t{gene_id}")

    if not records:
        raise SystemExit("No CDS sequences extracted - check the WARNINGs above")

    dest = outdir / args.out_name
    with dest.open("w") as fh:
        for header, seq in records:
            fh.write(header + "\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i : i + 70] + "\n")

    print(f"Wrote {dest} ({len(records)} transcript(s) from {n_genes_written} gene(s))")


if __name__ == "__main__":
    main()
