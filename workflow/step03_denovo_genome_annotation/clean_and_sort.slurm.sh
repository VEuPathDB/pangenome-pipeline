#!/usr/bin/env bash
#SBATCH --job-name=clean_sort
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }
# Step 03 works in one folder, $STEP03_OUT: the inputs you put there (genomes/genome.fna, the
# RNA-seq reads, curated_TE_library.fa) and every output. The paths below are relative to it.
mkdir -p "$STEP03_OUT" && cd "$STEP03_OUT"
# Environment & Shared paths (adjust to your cluster setup)
require_config FUNANNOTATE_SHARED
SHARED="$FUNANNOTATE_SHARED"   # shared funannotate install, set in pipeline_config.sh
GENOME_IN="genomes/genome.fna"
PREFIX="GENOME"
OUTDIR="funannotate"

mkdir -p "${OUTDIR}"

echo "Starting funannotate clean: $(date)"

singularity exec \
    -B "${PWD}:${PWD}" \
    -B "${SHARED}/overlay/usr/local/bin:/usr/local/bin" \
    -B "${SHARED}/overlay/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu" \
    --pwd "${PWD}" \
    --env PATH="/usr/local/bin:/venv/bin:/usr/bin:/bin" \
    --env GENEMARK_PATH="/usr/local/bin/genemark" \
    "${SHARED}/containers/funannotate.sif" \
    funannotate clean \
        -i "${GENOME_IN}" \
        -o "${OUTDIR}/${PREFIX}.clean.fna" \
        --minlen 500

echo "Starting funannotate sort: $(date)"

singularity exec \
    -B "${PWD}:${PWD}" \
    -B "${SHARED}/overlay/usr/local/bin:/usr/local/bin" \
    -B "${SHARED}/overlay/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu" \
    --pwd "${PWD}" \
    --env PATH="/usr/local/bin:/venv/bin:/usr/bin:/bin" \
    --env GENEMARK_PATH="/usr/local/bin/genemark" \
    "${SHARED}/containers/funannotate.sif" \
    funannotate sort \
        -i "${OUTDIR}/${PREFIX}.clean.fna" \
        -o "${OUTDIR}/${PREFIX}.sorted.fna" \
        -b "${PREFIX}"

echo "Finished clean and sort at $(date)"
