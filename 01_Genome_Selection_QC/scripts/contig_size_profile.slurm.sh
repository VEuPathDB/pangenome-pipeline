#!/bin/bash
#SBATCH --job-name=contig_size_profile
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=contig_size_profile_%j.out
#SBATCH --error=contig_size_profile_%j.err
#
# Per-genome contig size profile - how many contigs, how big, and how much
# of the assembly sits in short contigs (possible junk/contamination, not
# necessarily redundant with anything else - see contig_redundancy.py for
# the duplication screen, a different question).
#
# Output, all in OUTDIR: contig_sizes.tsv (every contig, + n_genes if GFF_INDIR is
# set), contig_size_summary.tsv (one row per genome), plots/contig_size_profile_by_genome.png
# (contig-count histogram by size bin, kb/Mb), plots/small_contig_count_by_genome.png, and
# (only if GFF_INDIR is set) plots/genes_by_contig_size_bin_by_genome.png,
# plots/pc_pseudogenes_by_contig_size_bin_by_genome.png, and small_contig_genes.tsv - every
# gene on a contig < MIN_LEN, with its GFF description, so mito/apicoplast genes (expected,
# not contamination) can be told apart from a real nuclear gene on a tiny contig (needs a
# closer look) by eye.
#
# Submit from your project folder (the one that holds genomes/):
#     sbatch 01_Genome_Selection_QC/scripts/contig_size_profile.slurm.sh
#
# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
GENOMES="genomes"                       # step-01 assemblies folder
GFF_INDIR="gff_originals"               # optional: one GFF3 per genome, for the gene-count plot.
                                         # Set to "" to skip gene counting (contig sizes only).
OUTDIR="contig_size_profile"            # where the tables/plots go
MIN_LEN=10000                           # contigs shorter than this (bp) count as "small"
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# Locate the .py. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there — fall back to the standard layout
# relative to the directory you submitted from.
PY=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/contig_size_profile.py" \
    "${SLURM_SUBMIT_DIR:-$PWD}/01_Genome_Selection_QC/scripts/contig_size_profile.py" \
    "01_Genome_Selection_QC/scripts/contig_size_profile.py"
do
    [[ -f "$cand" ]] && { PY="$cand"; break; }
done
[[ -n "$PY" ]] || { echo "ERROR: cannot find contig_size_profile.py (run from the project folder)" >&2; exit 1; }

[[ -d "$GENOMES" ]] || { echo "ERROR: genomes folder not found: $GENOMES (run from the project folder)" >&2; exit 1; }

### ADAPT ### — put samtools on the PATH (env stacking with the python env
#              below isn't reliable across micromamba versions, so .fai
#              generation and the python step use separate activations):
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate bioinfo
### END ADAPT ###

command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }

echo "=== 1/2: samtools faidx (skipping any genome that already has an up-to-date .fai) ==="
n=0
for ext in fna fasta fa; do
    for fasta in "$GENOMES"/*."$ext"; do
        [[ -e "$fasta" ]] || continue
        if [[ ! -f "$fasta.fai" || "$fasta" -nt "$fasta.fai" ]]; then
            samtools faidx "$fasta"
        fi
        n=$((n + 1))
    done
done
(( n > 0 )) || { echo "ERROR: no *.fna/*.fasta/*.fa files in $GENOMES" >&2; exit 1; }
echo "$n genome(s) indexed"

echo
echo "=== 2/2: contig_size_profile.py ==="
### ADAPT ### — put a python3 with pandas + matplotlib on the PATH:
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
python3 -c "import pandas, matplotlib" 2>/dev/null \
    || { echo "ERROR: python3 is missing pandas and/or matplotlib" >&2; exit 1; }

ARGS=(--genomes-dir "$GENOMES" --outdir "$OUTDIR" --min-len "$MIN_LEN")
[[ -n "${GFF_INDIR:-}" ]] && ARGS+=(--gff-indir "$GFF_INDIR")
python3 "$PY" "${ARGS[@]}"

echo
echo "Done — $OUTDIR/contig_sizes.tsv, $OUTDIR/contig_size_summary.tsv, $OUTDIR/plots/"
