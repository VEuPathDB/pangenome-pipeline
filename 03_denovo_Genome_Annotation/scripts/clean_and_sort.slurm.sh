#!/usr/bin/env bash
#SBATCH --job-name=clean_sort
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=clean_sort_%j.out
#SBATCH --error=clean_sort_%j.err

set -euo pipefail

# Environment & Shared paths (adjust to your cluster setup)
# Step 03 jobs are submitted from inside 03_denovo_Genome_Annotation/ (see its README), so look
# for the pipeline's lib/ upwards from there; PIPELINE_DIR overrides.
_pd="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: lib/slurm_common.sh not found above ${SLURM_SUBMIT_DIR:-$PWD} - set PIPELINE_DIR" >&2; exit 1; }
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
