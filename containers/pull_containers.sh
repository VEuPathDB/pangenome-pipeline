#!/usr/bin/env bash
#
# Download the tool images the pipeline's environments use, so jobs run without micromamba.
# Each environment `activate_env NAME` asks for becomes $CONTAINER_DIR/NAME.sif, or a folder
# $CONTAINER_DIR/NAME/ with one image per tool (lib/slurm_common.sh, activate_env).
#
#     bash containers/pull_containers.sh [CONTAINER_DIR]     # default: CONTAINER_DIR from pipeline_config.sh
#
# Images already there are skipped, so it can be re-run after adding one. Versions are pinned here;
# change a tag and delete the old .sif to upgrade. MCHelper has no published image - build it from
# containers/mchelper.def (see that file).

set -euo pipefail

source "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/../lib/slurm_common.sh"
CONTAINER_DIR="${1:-${CONTAINER_DIR:-}}"
[[ -n "$CONTAINER_DIR" ]] || _pipeline_die "give the folder: bash $0 /path/to/containers (or set CONTAINER_DIR in pipeline_config.sh)"

BIOCONTAINERS=docker://quay.io/biocontainers

# <file under CONTAINER_DIR>  <image>
IMAGES=(
    # python3 + pandas/matplotlib come from the pipeline's venv (uv sync --frozen --no-dev)
    "$ENV_SYNTENY/minimap2.sif" "$BIOCONTAINERS/minimap2:2.31--h118bc1c_0"
    "$ENV_SYNTENY/skani.sif" "$BIOCONTAINERS/skani:0.3.2--h79ce301_0"
    "$ENV_SYNTENY/samtools.sif" "$BIOCONTAINERS/samtools:1.24--h9dcdb79_1"
    "$ENV_SYNTENY/blast.sif" "$BIOCONTAINERS/blast:2.17.0--hb02a186_1"

    "$ENV_BIOINFO/agat.sif" "$BIOCONTAINERS/agat:1.7.0--pl5321hdfd78af_0"
    "$ENV_BIOINFO/samtools.sif" "$BIOCONTAINERS/samtools:1.24--h9dcdb79_1"
    "$ENV_BIOINFO/seqkit.sif" "$BIOCONTAINERS/seqkit:2.14.0--hb192632_0"

    "$ENV_REPEATMASKER.sif" "$BIOCONTAINERS/repeatmasker:4.2.4--pl5321hdfd78af_0"
    "$ENV_TE_ANNOT.sif" "$BIOCONTAINERS/te-aid:1.0.0--hdfd78af_0"
    "$ENV_GET_PANGENES.sif" "$BIOCONTAINERS/get_pangenes:20250904--hdfd78af_0"

    # step 02 TE library filtering (Step 2 / Step 4)
    "$ENV_TE_FILTER/blast.sif" "$BIOCONTAINERS/blast:2.17.0--hb02a186_1"
    "$ENV_TE_FILTER/cd-hit.sif" "$BIOCONTAINERS/cd-hit:4.8.1--hdbcaa40_0"
    "$ENV_TE_FILTER/seqkit.sif" "$BIOCONTAINERS/seqkit:2.14.0--hb192632_0"
)

if [[ -z "${CONTAINER_CMD:-}" ]]; then
    if command -v apptainer >/dev/null 2>&1; then CONTAINER_CMD=apptainer
    elif command -v singularity >/dev/null 2>&1; then CONTAINER_CMD=singularity
    else _pipeline_die "neither apptainer nor singularity found (try: module load apptainer)"
    fi
fi

failed=0
for ((i = 0; i < ${#IMAGES[@]}; i += 2)); do
    sif="$CONTAINER_DIR/${IMAGES[i]}"
    if [[ -s "$sif" ]]; then
        echo "have   ${IMAGES[i]}"
        continue
    fi
    mkdir -p "$(dirname "$sif")"
    echo "pull   ${IMAGES[i]}  <- ${IMAGES[i + 1]}"
    "$CONTAINER_CMD" pull --disable-cache "$sif" "${IMAGES[i + 1]}" || { echo "FAILED ${IMAGES[i]}" >&2; rm -f "$sif"; failed=1; }
done

[[ -s "$CONTAINER_DIR/$ENV_MCHELPER.sif" ]] \
    || echo "note   $ENV_MCHELPER.sif: build it with  $CONTAINER_CMD build --fakeroot $CONTAINER_DIR/$ENV_MCHELPER.sif $PIPELINE_DIR/containers/mchelper.def"
exit "$failed"
