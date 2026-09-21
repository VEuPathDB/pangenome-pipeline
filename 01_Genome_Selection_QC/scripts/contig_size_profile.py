#!/usr/bin/env python3
"""
Per-genome contig size profile: how many contigs, how big, and how much of
the assembly sits in short contigs that could be junk/contamination rather
than real sequence.

For each genome this script:
  1. Gets contig lengths from the assembly's .fai index (samtools faidx is
     run automatically if the index doesn't already exist next to the
     FASTA).
  2. Writes one summary row (n_contigs, total_bp, and how many
     contigs/how much bp fall below --min-len) and a full per-contig table.
  3. Plots each genome's contigs ranked largest-to-smallest, length on a
     log axis - the standard assembly-QC view for spotting the "cliff"
     between real chromosomes/scaffolds and a long tail of tiny contigs -
     plus one bar chart of how many small contigs each genome has.

This is a size/contamination-triage screen, not the redundancy screen
(contig_redundancy.py, which flags contigs duplicated elsewhere in the same
assembly) - a contig can be short and unique (this script) or short and a
near-duplicate of a longer contig (that one); they're different questions.

Usage:
    python3 contig_size_profile.py --genomes-dir genomes --outdir contig_size_profile
"""
import argparse
import collections
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

GENOME_EXTS = (".fna", ".fasta", ".fa")

# Same set/logic as contig_redundancy.py's gene_counts() - kept as its own
# local copy rather than importing across scripts, but must stay in sync:
# VEuPathDB/FungiDB-style GFFs use a specific feature type per biotype
# (protein_coding_gene, pseudogene, ...) instead of NCBI's generic "gene" +
# gene_biotype attribute, so both conventions need to be in this set.
GENE_FEATURE_TYPES = {
    "gene", "protein_coding_gene", "pseudogene", "ncRNA_gene",
    "tRNA_gene", "rRNA_gene", "snoRNA_gene", "snRNA_gene", "miRNA_gene",
}

# Same fallback logic as te_gene_overlap.py's infer_biotype(): an explicit
# gene_biotype/biotype/gene_type attribute wins when present (NCBI-style,
# "gene" in column 3 for everything); otherwise the feature type itself
# says the biotype (VEuPathDB/FungiDB-style, no biotype attribute at all).
GENE_FTYPE_TO_BIOTYPE = {
    "gene": "protein_coding",
    "protein_coding_gene": "protein_coding",
    "pseudogene": "pseudogene",
    "tRNA_gene": "tRNA",
    "rRNA_gene": "rRNA",
    "ncRNA_gene": "ncRNA",
    "snoRNA_gene": "snoRNA",
    "snRNA_gene": "snRNA",
    "miRNA_gene": "miRNA",
}
BIOTYPE_ATTR_KEYS = ("gene_biotype", "biotype", "gene_type")
PC_PSEUDO_BIOTYPES = {"protein_coding", "pseudogene"}  # what GET_PANGENES actually uses


def infer_biotype(ftype: str, attrs: dict) -> str:
    for key in BIOTYPE_ATTR_KEYS:
        if key in attrs:
            return attrs[key]
    return GENE_FTYPE_TO_BIOTYPE.get(ftype, "NA")


def parse_attributes(attr_str: str) -> dict:
    attrs = {}
    for kv in attr_str.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            attrs[k] = v
    return attrs


# Same order/keys as te_gene_overlap.py's gene_description(): NCBI puts
# "product" on the mRNA/CDS rather than the gene itself, so this alone
# won't catch every case - good enough for a quick "what is this gene"
# label without re-implementing the full mRNA-fallback logic here.
DESCRIPTION_ATTR_KEYS = ("description", "product", "Note", "note", "Name")


def gene_description(attrs: dict) -> str:
    from urllib.parse import unquote
    for key in DESCRIPTION_ATTR_KEYS:
        if key in attrs:
            return unquote(attrs[key]).replace("+", " ")
    return ""


def find_genomes(genomes_dir: Path) -> dict:
    """<genome_label>.fna/.fasta/.fa -> path, one per genome."""
    found = {}
    for ext in GENOME_EXTS:
        for p in sorted(genomes_dir.glob(f"*{ext}")):
            found.setdefault(p.stem, p)
    return found


def find_gff(gff_indir: Path, genome_label: str):
    """<genome_label>.gff/.gff3, falling back to a case-insensitive glob for
    naming-convention mismatches (mirrors the other scripts' find_gff)."""
    for ext in (".gff", ".gff3"):
        cand = gff_indir / f"{genome_label}{ext}"
        if cand.exists():
            return cand
    lowered = genome_label.lower()
    for f in gff_indir.glob("*.gff*"):
        if f.stem.lower() == lowered:
            return f
    return None


def parse_gene_records(gff_path: Path) -> list:
    """Every gene-level feature in a GFF3 as a list of dicts: contig,
    gene_id, ftype (the raw column-3 value - useful on its own as a
    mito/apicoplast tell, e.g. rRNA_gene clusters), biotype, description,
    start, end, strand. One pass over the file feeds both the per-contig
    gene counts and the small-contig gene-detail table, instead of parsing
    the GFF twice."""
    records = []
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] not in GENE_FEATURE_TYPES:
                continue
            attrs = parse_attributes(f[8])
            records.append({
                "contig": f[0],
                "gene_id": attrs.get("ID", ""),
                "ftype": f[2],
                "biotype": infer_biotype(f[2], attrs),
                "description": gene_description(attrs),
                "start": int(f[3]),
                "end": int(f[4]),
                "strand": f[6],
            })
    return records


def gene_counts_per_contig(records: list) -> tuple:
    """Returns (all_counts, pc_pseudo_counts) - both {seqid: n_genes} - from
    parse_gene_records()'s output. all_counts is every gene-level feature
    (same convention as contig_redundancy.py's gene_counts());
    pc_pseudo_counts is restricted to protein_coding/pseudogene - the
    biotypes GET_PANGENES actually uses, so a tRNA/rRNA/ncRNA gene sitting
    on an otherwise-empty small contig doesn't make that contig look
    load-bearing for the pangenome."""
    all_counts = collections.Counter()
    pc_pseudo_counts = collections.Counter()
    for rec in records:
        all_counts[rec["contig"]] += 1
        if rec["biotype"] in PC_PSEUDO_BIOTYPES:
            pc_pseudo_counts[rec["contig"]] += 1
    return all_counts, pc_pseudo_counts


def ensure_fai(fasta_path: Path) -> Path:
    """Returns the .fai path next to the FASTA, running `samtools faidx` if
    it doesn't already exist (or is older than the FASTA itself)."""
    fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
    if fai_path.exists() and fai_path.stat().st_mtime >= fasta_path.stat().st_mtime:
        return fai_path
    subprocess.run(["samtools", "faidx", str(fasta_path)], check=True)
    return fai_path


def read_fai(fai_path: Path) -> pd.DataFrame:
    """.fai columns: name, length, offset, linebases, linewidth."""
    df = pd.read_csv(fai_path, sep="\t", header=None,
                      names=["contig", "length", "offset", "linebases", "linewidth"])
    return df[["contig", "length"]]


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e6e6e0", zorder=0)
    ax.set_axisbelow(True)


SIZE_BINS_KB = [(0, 1), (1, 10), (10, 100), (100, 1000), (1000, float("inf"))]
SIZE_BIN_LABELS = ["<1 kb", "1-10 kb", "10-100 kb", "100 kb-1 Mb", ">=1 Mb"]
SIZE_BIN_COLORS = ["#C44E52", "#DD8452", "#CCB974", "#55A868", "#4C72B0"]


def size_bin_label(length_bp):
    kb = length_bp / 1000.0
    for (lo, hi), label in zip(SIZE_BINS_KB, SIZE_BIN_LABELS):
        if lo <= kb < hi:
            return label
    return SIZE_BIN_LABELS[-1]


def plot_size_histogram(per_contig: pd.DataFrame, dest: Path, ncols=3):
    """Small-multiples grid, one panel per genome: a histogram of contig
    count by size bin (kb, not bp - and categorical bins instead of a log
    axis, so there's nothing to "read" a log scale for). Bin edges
    (1 kb / 10 kb / 100 kb / 1 Mb) span the range a fungal/apicomplexan
    contig set actually takes - from stray short contigs up to whole
    chromosomes - while staying in human-readable kb/Mb units."""
    df = per_contig.copy()
    df["size_bin"] = df["length"].map(size_bin_label)
    genomes = sorted(df["genome_label"].unique())
    nrows = -(-len(genomes) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.0 * nrows), sharey=True)
    axes = axes.flatten()
    for ax, genome_label in zip(axes, genomes):
        sub = df[df["genome_label"] == genome_label]
        counts = sub["size_bin"].value_counts().reindex(SIZE_BIN_LABELS, fill_value=0)
        ax.bar(range(len(SIZE_BIN_LABELS)), counts.values, color=SIZE_BIN_COLORS, zorder=3)
        for x, v in enumerate(counts.values):
            if v > 0:
                ax.annotate(f"{int(v)}", xy=(x, v), xytext=(0, 2), textcoords="offset points",
                            ha="center", va="bottom", fontsize=7, color="#555")
        ax.set_xticks(range(len(SIZE_BIN_LABELS)))
        ax.set_xticklabels(SIZE_BIN_LABELS, rotation=45, ha="right", fontsize=7)
        style_axes(ax)
        ax.set_title(f"{genome_label} (n={len(sub)})", fontsize=10)
    for ax in axes[len(genomes):]:
        ax.axis("off")
    fig.supylabel("Number of contigs")
    fig.suptitle("Contig size distribution, by genome", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.96))
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_genes_by_size_bin(per_contig: pd.DataFrame, dest: Path, gene_col="n_genes",
                            label="genes", ncols=3):
    """Same layout as plot_size_histogram, but bars are total genes (from
    gene_col) on the contigs in that size bin, not contig count - directly
    answers "do the small contigs actually carry genes, or are they
    empty?" for each genome. Only meaningful if gene_col was populated
    (--gff-indir given). gene_col="n_genes" = every gene-level feature;
    "n_pc_pseudo" = protein_coding/pseudogene only, i.e. what GET_PANGENES
    actually uses - a tRNA sitting alone on a small contig makes that
    contig "load-bearing" in n_genes but not in n_pc_pseudo."""
    df = per_contig.copy()
    df["size_bin"] = df["length"].map(size_bin_label)
    genomes = sorted(df["genome_label"].unique())
    nrows = -(-len(genomes) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols, 3.0 * nrows), sharey=True)
    axes = axes.flatten()
    for ax, genome_label in zip(axes, genomes):
        sub = df[df["genome_label"] == genome_label]
        gene_sums = sub.groupby("size_bin")[gene_col].sum().reindex(SIZE_BIN_LABELS, fill_value=0)
        ax.bar(range(len(SIZE_BIN_LABELS)), gene_sums.values, color=SIZE_BIN_COLORS, zorder=3)
        for x, v in enumerate(gene_sums.values):
            if v > 0:
                ax.annotate(f"{int(v)}", xy=(x, v), xytext=(0, 2), textcoords="offset points",
                            ha="center", va="bottom", fontsize=7, color="#555")
        ax.set_xticks(range(len(SIZE_BIN_LABELS)))
        ax.set_xticklabels(SIZE_BIN_LABELS, rotation=45, ha="right", fontsize=7)
        style_axes(ax)
        ax.set_title(f"{genome_label} (n={int(sub[gene_col].sum())} {label})", fontsize=10)
    for ax in axes[len(genomes):]:
        ax.axis("off")
    fig.supylabel(f"Number of {label}")
    fig.suptitle(f"Annotated {label} by contig size bin, by genome", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.96))
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_small_contig_count(summary: pd.DataFrame, min_len: int, dest: Path):
    """One horizontal bar per genome: how many contigs fall below --min-len."""
    values = summary.set_index("genome_label")["n_contigs_below"].sort_values()
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(values) + 1.5))
    ax.barh(values.index, values.values, color="#DD8452", height=0.65, zorder=3)
    style_axes(ax)
    for y, v in enumerate(values.values):
        ax.annotate(f"{int(v)}", xy=(v, y), xytext=(5, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=8, color="#555")
    ax.set_xlabel(f"Contigs < {min_len:,} bp")
    ax.set_title(f"Contigs under {min_len:,} bp, by genome", fontsize=13, fontweight="bold", pad=12)
    ax.text(0.98, 0.03, f"total = {int(values.sum()):,}", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9, color="#555")
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes-dir", required=True,
                     help="Directory with one assembly FASTA per genome (<genome_label>.fna/.fasta/.fa)")
    ap.add_argument("--outdir", required=True, help="Directory to write the summary/plots to")
    ap.add_argument("--min-len", type=int, default=10000,
                     help="Contigs shorter than this (bp) count toward n_contigs_below/bp_below (default: 10000)")
    ap.add_argument("--gff-indir", default=None,
                     help="Optional: directory with one GFF3 per genome (<genome_label>.gff/.gff3). "
                          "If given, adds a gene count per contig and a genes-by-size-bin plot.")
    args = ap.parse_args()

    genomes_dir = Path(args.genomes_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    gff_indir = Path(args.gff_indir) if args.gff_indir else None

    genomes = find_genomes(genomes_dir)
    if not genomes:
        raise SystemExit(f"No {'/'.join(GENOME_EXTS)} files found in {genomes_dir}")

    if subprocess.run(["which", "samtools"], capture_output=True).returncode != 0:
        print("WARNING: samtools not found on PATH - this only works if every "
              "genome already has an up-to-date .fai next to it", file=sys.stderr)

    all_rows = []
    summary_rows = []
    small_contig_gene_rows = []
    for genome_label, fasta_path in sorted(genomes.items()):
        try:
            fai_path = ensure_fai(fasta_path)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"WARNING: could not get/build .fai for {genome_label} ({fasta_path.name}): {e}",
                  file=sys.stderr)
            continue
        contigs = read_fai(fai_path)
        contigs.insert(0, "genome_label", genome_label)

        gene_records = []
        if gff_indir is not None:
            gff_path = find_gff(gff_indir, genome_label)
            if gff_path is None:
                print(f"WARNING: no GFF found for {genome_label} in {gff_indir} - "
                      f"n_genes will be 0 for all its contigs", file=sys.stderr)
                contigs["n_genes"] = 0
                contigs["n_pc_pseudo"] = 0
            else:
                gene_records = parse_gene_records(gff_path)
                all_counts, pc_pseudo_counts = gene_counts_per_contig(gene_records)
                contigs["n_genes"] = contigs["contig"].map(all_counts).fillna(0).astype(int)
                contigs["n_pc_pseudo"] = contigs["contig"].map(pc_pseudo_counts).fillna(0).astype(int)
        else:
            contigs["n_genes"] = 0
            contigs["n_pc_pseudo"] = 0
        all_rows.append(contigs)

        below = contigs[contigs["length"] < args.min_len]

        if gene_records:
            small_contig_lengths = dict(zip(below["contig"], below["length"]))
            for rec in gene_records:
                contig_len = small_contig_lengths.get(rec["contig"])
                if contig_len is None:
                    continue  # gene's contig isn't in the small-contig set
                small_contig_gene_rows.append({
                    "genome_label": genome_label,
                    "contig": rec["contig"],
                    "contig_length_bp": contig_len,
                    "gene_id": rec["gene_id"],
                    "feature_type": rec["ftype"],
                    "biotype": rec["biotype"],
                    "description": rec["description"],
                    "start": rec["start"],
                    "end": rec["end"],
                    "strand": rec["strand"],
                })
        total_bp = int(contigs["length"].sum())
        bp_below = int(below["length"].sum())
        n_genes_below = int(below["n_genes"].sum())
        n_pc_pseudo_below = int(below["n_pc_pseudo"].sum())
        summary_rows.append({
            "genome_label": genome_label,
            "n_contigs": len(contigs),
            "total_bp": total_bp,
            "n_contigs_below": len(below),
            "bp_below": bp_below,
            "pct_bp_below": round(100.0 * bp_below / total_bp, 3) if total_bp else 0.0,
            "n_genes_below": n_genes_below,
            "n_pc_pseudo_below": n_pc_pseudo_below,
            "largest_contig_bp": int(contigs["length"].max()) if len(contigs) else 0,
            "smallest_contig_bp": int(contigs["length"].min()) if len(contigs) else 0,
        })
        print(f"{genome_label}: {len(contigs)} contigs, {len(below)} below {args.min_len:,} bp "
              f"({bp_below:,} / {total_bp:,} bp = {summary_rows[-1]['pct_bp_below']}%), "
              f"{n_genes_below} gene(s) ({n_pc_pseudo_below} protein_coding/pseudogene) on those small contigs")

    if not all_rows:
        raise SystemExit("No genome could be indexed - nothing to summarize")

    per_contig = pd.concat(all_rows, ignore_index=True)
    per_contig.to_csv(outdir / "contig_sizes.tsv", sep="\t", index=False)

    summary = pd.DataFrame(summary_rows).sort_values("pct_bp_below", ascending=False)
    summary.to_csv(outdir / "contig_size_summary.tsv", sep="\t", index=False)
    print(f"\nWrote {outdir / 'contig_sizes.tsv'} ({len(per_contig)} contigs total) and "
          f"{outdir / 'contig_size_summary.tsv'} ({len(summary)} genomes)")

    if gff_indir is not None:
        small_genes = pd.DataFrame(small_contig_gene_rows)
        small_genes.to_csv(outdir / "small_contig_genes.tsv", sep="\t", index=False)
        print(f"Wrote {outdir / 'small_contig_genes.tsv'} ({len(small_genes)} gene(s) on "
              f"contigs < {args.min_len:,} bp - check these by hand: mito/apicoplast genes "
              f"are expected here and aren't contamination, a real nuclear gene on a tiny "
              f"contig needs a closer look, anything that reads as foreign/host DNA is the "
              f"actual contamination signal)")

    plots_dir = outdir / "plots"
    plots_dir.mkdir(exist_ok=True)
    plot_size_histogram(per_contig, plots_dir / "contig_size_profile_by_genome.png")
    plot_small_contig_count(summary, args.min_len, plots_dir / "small_contig_count_by_genome.png")
    if gff_indir is not None:
        plot_genes_by_size_bin(per_contig, plots_dir / "genes_by_contig_size_bin_by_genome.png",
                                gene_col="n_genes", label="genes")
        plot_genes_by_size_bin(per_contig, plots_dir / "pc_pseudogenes_by_contig_size_bin_by_genome.png",
                                gene_col="n_pc_pseudo", label="protein_coding/pseudogene genes")
    print(f"Plots written to {plots_dir}/")


if __name__ == "__main__":
    main()
