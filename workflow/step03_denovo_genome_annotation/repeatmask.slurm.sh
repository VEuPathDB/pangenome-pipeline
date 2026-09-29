#!/usr/bin/env bash
#SBATCH --job-name=repeatmask
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
require_config FUNANNOTATE_SHARED
SHARED="$FUNANNOTATE_SHARED"   # shared funannotate install, set in pipeline_config.sh
GENOME_IN="funannotate/GENOME.sorted.fna"
TE_LIB="curated_TE_library.fa"
OUTDIR="funannotate/repeatmask"

mkdir -p "${OUTDIR}"

echo "Starting RepeatMasker: $(date)"

# Option A: Using RepeatMasker directly
singularity exec \
    -B "${PWD}:${PWD}" \
    "${SHARED}/containers/dfam-tetools-latest.sif" \
    RepeatMasker \
        -pa 32 \
        -lib "${TE_LIB}" \
        -dir "${OUTDIR}" \
        -xsmall \
        -gff \
        -html \
        "${GENOME_IN}"

# Option B (alternative): Using funannotate mask
# singularity exec \
#     -B "${PWD}:${PWD}" \
#     -B "${SHARED}/databases/funannotate_db:${SHARED}/databases/funannotate_db" \
#     --env FUNANNOTATE_DB="${SHARED}/databases/funannotate_db" \
#     "${SHARED}/containers/funannotate.sif" \
#     funannotate mask \
#         -i "${GENOME_IN}" \
#         -o "${OUTDIR}/GENOME.masked.fna" \
#         --cpus 32

echo "Finished repeat masking at $(date)"
