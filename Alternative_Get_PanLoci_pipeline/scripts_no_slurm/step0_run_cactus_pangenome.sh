#!/bin/bash
set -euo pipefail

# ==============================================================================
# STEP 0: CACTUS PANGENOME RUNNER (NON-SLURM)
# ==============================================================================

# Main directory where Cactus will run.
PROJECT_DIR="/path/to/panlocus_run"

# Cactus input genome list.
# Format: one genome per line, usually:
#   strain_id<TAB>/absolute/path/to/genome.fasta
GENOMES_FILE="genomes_all.txt"

# Cactus Apptainer/Singularity image.
CACTUS_SIF="/path/to/cactus.sif"

# Cactus output configuration.
OUT_DIR="cactus_output"
OUT_NAME="pangenome"
JOBSTORE="jobStore"

# Cactus reference genome name.
CACTUS_REFERENCE="REF3D7"

# Resources passed to cactus-pangenome.
MAX_CORES="8"
MAX_MEMORY="16G"

# Bind mount needed by Apptainer/Singularity so absolute storage paths resolve
# correctly inside the container.
APPTAINER_BIND_PATH="/path/to/storage"

# Extra Cactus options used for this analysis.
CACTUS_EXTRA_OPTIONS="--vcf --giraffe --gfa --gbz --viz --permissiveContigFilter"

# ==============================================================================
# ENVIRONMENT AND INPUT CHECKS
# ==============================================================================

cd "${PROJECT_DIR}"

if [[ ! -f "${GENOMES_FILE}" ]]; then
    echo "[ERROR] Genome list not found: ${PROJECT_DIR}/${GENOMES_FILE}" >&2
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
# Running Cactus pangenome alignment
echo "Project directory:     ${PROJECT_DIR}"
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
echo "${PROJECT_DIR}/${OUT_DIR}/${OUT_NAME}.full.hal"
