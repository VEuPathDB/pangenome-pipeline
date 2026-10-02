#!/bin/bash
#SBATCH --job-name=RepeatMasker
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Annotate TEs in every assembly with RepeatMasker, using the curated TE library.
# One output folder per genome: $STEP02_OUT/repeatmasker/RM_<genome_id>/  (genome_id = the
# .fna file name without the extension; an NCBI accession is trimmed to GCA_/GCF_).
#
# Submit from your project folder (the one that holds genomes/):
#     sbatch $PIPELINE_DIR/workflow/step02_te_annotation_clean_up/run_repeatmasker.slurm.sh
#
set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
GENOMES="genomes"                     # <genome_id>.fna assemblies (step-01 folder), or a full path
LIB="$STEP02_OUT/curated_TE_library.fa"   # curated TE library FASTA (put it in the step-02 output folder), or a full path
RM_OUTDIR="$STEP02_OUT/repeatmasker"  # one RM_<genome_id>/ per genome goes here
# ────────────────────────────────────────────────────────────────────────────

### ADAPT ### — put RepeatMasker on the PATH (must be configured with a Dfam
#              library; see the README if `famdb.py … FamDB data directory not found`):
activate_env "$ENV_REPEATMASKER"  # env names / backend: pipeline_config.sh
### END ADAPT ###

command -v RepeatMasker >/dev/null || { echo "ERROR: RepeatMasker not on PATH" >&2; exit 1; }
[[ -d "$GENOMES" ]] || { echo "ERROR: genomes folder not found: $GENOMES (run from the project folder)" >&2; exit 1; }
[[ -s "$LIB" ]]     || { echo "ERROR: TE library not found: $LIB" >&2; exit 1; }

CPUS="${SLURM_CPUS_PER_TASK:-24}"
PA=$(( CPUS / 4 )); (( PA > 0 )) || PA=1        # rmblast uses 4 cores per -pa unit

for fna in "$GENOMES"/*.fna; do
    [[ -e "$fna" ]] || { echo "ERROR: no *.fna in $GENOMES" >&2; exit 1; }

    name=$(basename "$fna" .fna)
    acc=$(echo "$name" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' || true)
    [[ -n "$acc" ]] && name="$acc"

    outdir="$RM_OUTDIR/RM_${name}"
    mkdir -p "$outdir"

    echo "=========================================================="
    echo "Processing: $(basename "$fna")  ->  $outdir/   (-pa $PA)"
    echo "=========================================================="

    RepeatMasker \
        -pa "$PA" \
        -lib "$LIB" \
        -dir "$outdir" \
        -no_is \
        -gff \
        -xsmall \
        -html \
        "$fna"

    echo "Finished $name at $(date)"
done

echo "All RepeatMasker jobs completed successfully!"
