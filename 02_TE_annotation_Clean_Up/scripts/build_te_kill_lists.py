#!/usr/bin/env python3
"""
Build a single AGAT "kill list" of TE-derived gene IDs, to strip those genes
back out of every Liftoff-lifted GFF that carries them (Liftoff propagates a
source annotation's gene calls - including bad ones - to every genome it
gets lifted onto).

Works straight off gene_te_overlap.tsv (te_gene_overlap.py's output) - no
need to have run reclassify_te_genes.py first. Selects every protein_coding
gene with TE_overlap_pct >= --min-pct (default: 80, the same "likely a TE"
cutoff used elsewhere in this project - override it here if needed without
touching any other script/output).

One combined file, not one per genome: a raw GFF ID like "gene-MGG_02216" is
only ever meaningful within its own genome's GFF, and a single Liftoff-
lifted file only ever contains genes lifted from one source annotation (the
normal Liftoff output shape) - so within any one lifted file, an ID can only
match the one gene it's actually supposed to. That also means you can run
the exact same kill list against every lifted GFF, whichever source genome
it came from: AGAT simply won't find IDs that aren't present in a given
file, so nothing is removed from a file that doesn't carry any flagged
genes. (This stops being safe only if lifted output ever mixes genes from
several different source annotations into one GFF file - not how Liftoff
normally writes it.)

This script never opens a GFF (works purely off gene_te_overlap.tsv), so it
has no way to know which ID convention each gene's *source* annotation
actually used - NCBI-style GFFs prefix gene IDs with "gene-" (which
te_gene_overlap.py's strip_prefix removes to get the plain gene_id stored in
the table), VEuPathDB/FungiDB GFFs (and Liftoff output lifted from one)
generally don't prefix at all, and Ensembl-style ones sometimes use "gene:".
Liftoff carries the source annotation's ID straight through, so guessing
wrong here would silently leave every TE-derived gene lifted from a
non-NCBI-style source in place. Rather than guess, all three forms are
written for every gene_id; the extra lines that don't match anything in a
given lifted GFF are simply no-ops for AGAT (see above).

Usage:
    python3 build_te_kill_lists.py \
        --gene-table te_annotation_integration/gene_te_overlap.tsv \
        --outdir kill_lists \
        [--min-pct 80] [--extra-genes confirmed_20_80.tsv]

Writes <outdir>/te_derived_genes.kill_list.txt (one "gene-<gene_id>" per
line), ready for:

    agat_sp_filter_feature_from_kill_list.pl \
        --gff <lifted.gff> \
        --kill_list kill_lists/te_derived_genes.kill_list.txt \
        --output <lifted.filtered.gff>

AGAT matches the "ID" attribute by default and removing a gene
automatically cascades to its mRNA/exon/CDS children, so no extra flags are
needed.

--extra-genes takes the same 2-column TSV (genome_label<TAB>gene_id, no
header) used by reclassify_te_genes.py/extract_flagged_genes.py, for genes
confirmed manually from the 20-80% inspect band that should be killed too.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd


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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-table", required=True,
                     help="gene_te_overlap.tsv produced by te_gene_overlap.py")
    ap.add_argument("--outdir", required=True,
                     help="Directory to write te_derived_genes.kill_list.txt to")
    ap.add_argument("--min-pct", type=float, default=80.0,
                     help="Minimum TE_overlap_pct to include a protein_coding gene (default: 80)")
    ap.add_argument("--extra-genes", default=None,
                     help="Optional 2-column TSV (genome_label<TAB>gene_id, no header) of "
                          "additional genes to include on top of --min-pct")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.gene_table, sep="\t")
    auto = df[(df["gene_category"] == "protein_coding") & (df["TE_overlap_pct"] >= args.min_pct)]
    pairs = set(zip(auto["genome_label"], auto["gene_id"]))
    print(f"{len(pairs)} gene(s) with TE_overlap_pct >= {args.min_pct}")

    if args.extra_genes:
        extra_pairs = load_extra_genes(args.extra_genes)
        valid_pairs = set(zip(df["genome_label"], df["gene_id"]))
        n_added = 0
        for pair in extra_pairs:
            if pair not in valid_pairs:
                print(f"WARNING: --extra-genes entry {pair} not found in {args.gene_table}, skipping")
                continue
            if pair not in pairs:
                n_added += 1
            pairs.add(pair)
        print(f"Added {n_added} gene(s) from {args.extra_genes}")

    if not pairs:
        raise SystemExit(f"No genes selected (TE_overlap_pct >= {args.min_pct}, "
                          f"plus none/no valid entries in --extra-genes)")

    # Same raw ID showing up under two different genomes would be ambiguous
    # if a lifted file ever mixed genes from multiple sources - flag it so
    # it's not a silent surprise, even though it's harmless for normal
    # one-source-per-file Liftoff output.
    genomes_by_gene_id = defaultdict(set)
    for genome_label, gene_id in pairs:
        genomes_by_gene_id[gene_id].add(genome_label)
    dupes = {gid: genomes for gid, genomes in genomes_by_gene_id.items() if len(genomes) > 1}
    if dupes:
        print(f"NOTE: {len(dupes)} gene ID(s) are shared by more than one genome_label "
              f"(e.g. {next(iter(dupes.items()))}) - fine as long as no lifted GFF mixes "
              f"genes lifted from multiple source annotations into one file.")

    gene_ids = sorted({gene_id for _, gene_id in pairs})
    dest = outdir / "te_derived_genes.kill_list.txt"
    # Write every ID-prefix convention seen in this project (bare, NCBI's
    # "gene-", Ensembl-style "gene:") for each gene - see the module
    # docstring for why guessing a single one isn't safe. AGAT ignores kill
    # list entries that don't match anything in a given GFF.
    with dest.open("w") as fh:
        for gene_id in gene_ids:
            for candidate in (gene_id, f"gene-{gene_id}", f"gene:{gene_id}"):
                fh.write(candidate + "\n")

    print(f"\n{len(gene_ids)} unique gene ID(s) from {len(pairs)} (genome, gene) pair(s) "
          f"written to {dest} (3 ID-convention variants per gene)")


if __name__ == "__main__":
    main()
