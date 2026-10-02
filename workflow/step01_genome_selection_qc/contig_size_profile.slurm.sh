#!/bin/bash
#SBATCH --job-name=contig_size_profile
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Per-genome contig size profile - how many contigs, how big, and how much
# of the assembly sits in short contigs (possible junk/contamination, not
# necessarily redundant with anything else - see contig_redundancy.py for
# the duplication screen, a different question).
#
# Output, all in OUTDIR ($STEP01_OUT/contig_size_profile): contig_sizes.tsv (every contig, + n_genes if GFF_INDIR is
# set), contig_size_summary.tsv (one row per genome), plots/contig_size_profile_by_genome.png
# (contig-count histogram by size bin, kb/Mb), plots/small_contig_count_by_genome.png, and
# (only if GFF_INDIR is set) plots/genes_by_contig_size_bin_by_genome.png,
# plots/pc_pseudogenes_by_contig_size_bin_by_genome.png, and small_contig_genes.tsv - every
# gene on a contig < MIN_LEN, with its GFF description, so mito/apicoplast genes (expected,
# not contamination) can be told apart from a real nuclear gene on a tiny contig (needs a
# closer look) by eye.
#
# Submit from your project folder (the one that holds genomes/):
#     sbatch $PIPELINE_DIR/workflow/step01_genome_selection_qc/contig_size_profile.slurm.sh
#
set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
GENOMES="genomes"                       # step-01 assemblies folder
GFF_INDIR="gff_originals"               # optional: one GFF3 per genome, for the gene-count plot.
                                         # Set to "" to skip gene counting (contig sizes only).
OUTDIR="$STEP01_OUT/contig_size_profile"   # where the tables/plots go
MIN_LEN=10000                           # contigs shorter than this (bp) count as "small"
# ────────────────────────────────────────────────────────────────────────────

[[ -d "$GENOMES" ]] || { echo "ERROR: genomes folder not found: $GENOMES (run from the project folder)" >&2; exit 1; }

### ADAPT ### — put samtools on the PATH (env stacking with the python env
#              below isn't reliable across micromamba versions, so .fai
#              generation and the python step use separate activations):
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
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
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
python3 -c "import pandas, matplotlib" 2>/dev/null \
    || { echo "ERROR: python3 is missing pandas and/or matplotlib" >&2; exit 1; }

ARGS=(--genomes-dir "$GENOMES" --outdir "$OUTDIR" --min-len "$MIN_LEN")
[[ -n "${GFF_INDIR:-}" ]] && ARGS+=(--gff-indir "$GFF_INDIR")
python3 -m pangenome_pipeline.step01_genome_selection_qc.contig_size_profile "${ARGS[@]}"

echo
echo "Done — $OUTDIR/contig_sizes.tsv, $OUTDIR/contig_size_summary.tsv, $OUTDIR/plots/"
