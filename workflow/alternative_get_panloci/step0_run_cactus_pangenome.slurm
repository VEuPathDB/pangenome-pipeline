#!/bin/bash
#SBATCH --job-name=cactus_pangenome
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=200:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=low

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT (ALT_PANLOCI_OUT), pipeline_config.sh, src/ on PYTHONPATH
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ==============================================================================
# STEP 0: CACTUS PANGENOME CONFIGURATION
# ==============================================================================

# Main directory where Cactus will run. This directory should contain the
# genomes file listed below.
# Default: $ALT_PANLOCI_OUT (output/alternative_get_panloci/ in the project folder), the same
# folder run_get_panloci uses; set RUN_DIR to use another.
RUN_DIR="${RUN_DIR:-$ALT_PANLOCI_OUT}"

# Cactus input genome list.
# Format: one genome per line, usually:
#   strain_id<TAB>/absolute/path/to/genome.fasta
# The strain_id values must match the names later used in samples.tsv.
GENOMES_FILE="genomes_all.txt"

# Cactus Apptainer/Singularity image.
CACTUS_SIF="/path/to/cactus.sif"

# Cactus output configuration.
OUT_DIR="cactus_output"
OUT_NAME="pangenome"
JOBSTORE="jobStore"

# Cactus reference genome name.
# This is the Cactus graph/reference choice. It does not have to be the same as
# the panlocus seed strain, which is controlled later by samples.tsv.
CACTUS_REFERENCE="REF3D7"

# Resources passed to cactus-pangenome.
MAX_CORES="32"
MAX_MEMORY="60G"

# Bind mount needed by Apptainer/Singularity so absolute storage paths resolve
# correctly inside the container.
APPTAINER_BIND_PATH="/path/to/storage"

# Extra Cactus options used for this analysis.
# --permissiveContigFilter allows permissive handling of shorter contigs.
CACTUS_EXTRA_OPTIONS="--vcf --giraffe --gfa --gbz --viz --permissiveContigFilter"

# ==============================================================================
# ENVIRONMENT AND INPUT CHECKS
# ==============================================================================

mkdir -p "${RUN_DIR}"
cd "${RUN_DIR}"

if [[ ! -f "${GENOMES_FILE}" ]]; then
    echo "[ERROR] Genome list not found: ${RUN_DIR}/${GENOMES_FILE}" >&2
    exit 1
fi

if [[ ! -f "${CACTUS_SIF}" ]]; then
    echo "[ERROR] Cactus SIF image not found: ${CACTUS_SIF}" >&2
    exit 1
fi

mkdir -p ./apptainer_tmp ./apptainer_cache "${OUT_DIR}"
export APPTAINER_TMPDIR="$(pwd)/apptainer_tmp"
export APPTAINER_CACHEDIR="$(pwd)/apptainer_cache"
export APPTAINER_BIND="${APPTAINER_BIND_PATH}"

echo "======================================================================"
echo " Step 0: Cactus pangenome alignment"
echo "======================================================================"
echo "Run directory:         ${RUN_DIR}"
echo "Genome list:           ${GENOMES_FILE}"
echo "Cactus image:          ${CACTUS_SIF}"
echo "Output directory:      ${OUT_DIR}"
echo "Output name:           ${OUT_NAME}"
echo "JobStore:              ${JOBSTORE}"
echo "Cactus reference:      ${CACTUS_REFERENCE}"
echo "Max cores:             ${MAX_CORES}"
echo "Max memory:            ${MAX_MEMORY}"
echo "APPTAINER_BIND:        ${APPTAINER_BIND}"
echo "Extra options:         ${CACTUS_EXTRA_OPTIONS}"
echo "======================================================================"

# ==============================================================================
# RUN CACTUS
# ==============================================================================

apptainer exec "${CACTUS_SIF}" \
    cactus-pangenome "./${JOBSTORE}" "${GENOMES_FILE}" \
    --outDir "./${OUT_DIR}" \
    --outName "${OUT_NAME}" \
    --reference "${CACTUS_REFERENCE}" \
    ${CACTUS_EXTRA_OPTIONS} \
    --maxCores "${MAX_CORES}" \
    --maxMemory "${MAX_MEMORY}" \
    --workDir "./apptainer_tmp"

echo "======================================================================"
echo " Step 0 completed"
echo "======================================================================"
echo "Main HAL output:"
echo "${RUN_DIR}/${OUT_DIR}/${OUT_NAME}.full.hal"
