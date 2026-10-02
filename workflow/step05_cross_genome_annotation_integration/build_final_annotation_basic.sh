#!/bin/bash
#SBATCH --job-name=build_final_annotation
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G                   # AGAT holds every file in memory; round 1 combines all other genomes at once
#SBATCH --time=48:00:00

# Submit from your project folder (logs/ must exist there - SLURM does not create it, and
# check_final_annotations_job.sh reads these jobs' logs from it). Either every genome one after
# another in a single job:
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/build_final_annotation_basic.sh
# or one job per genome, in parallel (N = number of genomes in ORDER):
#     N=$(grep -c . genome_order.txt)
#     sbatch --array=1-${N} $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/build_final_annotation_basic.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
INPUT_DIR="$STEP05_OUT/gff_for_liftoff"   # each genome's own annotation, used to fill what Liftoff did not
                                            # transfer (incl. the non-coding genes) - the prepare_for_liftoff output
PAIRWISE_DIR="$STEP05_OUT/liftoff_output/polished"   # directory with your all-against-all Liftoff output
FINAL_DIR="$STEP05_OUT/expanded_annotation"   # where to write this pipeline's output
ORDER="genome_order.txt"                   # first genome = reference; comma-separated list or a
                                            # path to a file with one genome name per line
ORIGINAL_SUFFIX=".gff3"                    # filename suffix for each genome's own GFF in INPUT_DIR
PAIR_SUFFIX="_liftoff.gff3_polished"       # filename suffix your Liftoff output has
MIN_COVERAGE="0.9"                         # minimum Liftoff coverage for a "good" (round-1) transfer
MIN_IDENTITY="0.9"                         # minimum Liftoff sequence identity for a "good" (round-1) transfer

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (change the env name if your
### cluster's AGAT environment isn't called this):
activate_env "$ENV_BIOINFO"  # env names / backend: pipeline_config.sh
### END ADAPT ###
set -u

find_script BUILD_SCRIPT workflow/step05_cross_genome_annotation_integration/build_final_annotation.sh

# As a SLURM array, task N builds only the Nth genome of ORDER.
TARGET_ARGS=()
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    order_list=()
    if [[ -f "$ORDER" ]]; then
        while IFS= read -r line; do
            line="${line%%#*}"; line="${line//[[:space:]]/}"
            [[ -n "$line" ]] && order_list+=("$line")
        done < "$ORDER"
    else
        IFS=',' read -ra order_list <<< "$ORDER"
    fi
    idx=$((SLURM_ARRAY_TASK_ID - 1))
    [[ $idx -ge 0 && $idx -lt ${#order_list[@]} ]] || { echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} is outside 1-${#order_list[@]}" >&2; exit 1; }
    TARGET_ARGS=(-t "${order_list[$idx]}")
    echo "Array task ${SLURM_ARRAY_TASK_ID}: building ${order_list[$idx]}"
fi

"$BUILD_SCRIPT" \
  -i "$INPUT_DIR" -p "$PAIRWISE_DIR" -o "$FINAL_DIR" -r "$ORDER" \
  -g "$ORIGINAL_SUFFIX" -x "$PAIR_SUFFIX" \
  -c "$MIN_COVERAGE" -d "$MIN_IDENTITY" \
  ${TARGET_ARGS[@]+"${TARGET_ARGS[@]}"}

echo "Build final annotation completed. Outputs are in $FINAL_DIR"
