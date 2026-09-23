#!/bin/bash
#SBATCH --job-name=build_final_annotation
#SBATCH --partition=low
#SBATCH --output=expanded_annotation/logs/%x_%j.out
#SBATCH --error=expanded_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G                   # AGAT holds every file in memory; round 1 combines all other genomes at once
#SBATCH --time=48:00:00

# Submit from your project folder. The log folder must exist BEFORE submitting (SLURM does not
# create it - without it the job leaves no log), and it must match the --output lines above.
# Either every genome one after another in a single job:
#     mkdir -p expanded_annotation/logs && sbatch 05_Cross-genome_annotation_integration/scripts/build_final_annotation_basic.sh
# or one job per genome, in parallel (N = number of genomes in ORDER):
#     mkdir -p expanded_annotation/logs
#     N=$(grep -c . genome_order.txt)
#     sbatch --array=1-${N} 05_Cross-genome_annotation_integration/scripts/build_final_annotation_basic.sh

set -euo pipefail

# Edit these for your run.
INPUT_DIR="gff_for_liftoff"                # each genome's own annotation, used to fill what Liftoff did not
                                            # transfer (incl. the non-coding genes) - the prepare_for_liftoff output
PAIRWISE_DIR="liftoff_output/polished"     # directory with your all-against-all Liftoff output
OUTPUT_DIR="expanded_annotation"           # where to write this pipeline's output (logs go to expanded_annotation/logs)
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
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate bioinfo
### END ADAPT ###
set -u

# Locate the script. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there - fall back to the standard layout
# relative to the directory you submitted from.
BUILD_SCRIPT=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/build_final_annotation.sh" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Cross-genome_annotation_integration/scripts/build_final_annotation.sh" \
    "05_Cross-genome_annotation_integration/scripts/build_final_annotation.sh"
do
    [[ -f "$cand" ]] && { BUILD_SCRIPT="$cand"; break; }
done
[[ -n "$BUILD_SCRIPT" ]] || { echo "ERROR: cannot find build_final_annotation.sh (run from the project folder)" >&2; exit 1; }

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
  -i "$INPUT_DIR" -p "$PAIRWISE_DIR" -o "$OUTPUT_DIR" -r "$ORDER" \
  -g "$ORIGINAL_SUFFIX" -x "$PAIR_SUFFIX" \
  -c "$MIN_COVERAGE" -d "$MIN_IDENTITY" \
  ${TARGET_ARGS[@]+"${TARGET_ARGS[@]}"}

echo "Build final annotation completed. Outputs are in $OUTPUT_DIR"
