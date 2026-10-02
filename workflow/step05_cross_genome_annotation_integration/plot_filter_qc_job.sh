#!/bin/bash
#SBATCH --job-name=plot_filter_qc
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

# Run this after 7a (filter_final_annotation_job.sh) has finished. Submit from your project
# folder (logs/ must exist there - SLURM does not create it):
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/plot_filter_qc_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
STRAINS="genome_order.txt"                   # comma-separated names or a file with one per line, first = the reference (same as ORDER)
ORIGINAL_DIR="$STEP05_OUT/expanded_annotation"   # build_final_annotation_basic.sh's FINAL_DIR (before filtering)
FILTERED_DIR="$STEP05_OUT/filtered_annotation"   # 7a's OUT_DIR (after filtering)
MIN_OVERLAP_BP="20"                          # same --min-overlap-bp used for the real 7a run
OUT_DIR="$STEP05_OUT/filter_qc_plots"   # where to write the plots and summary.tsv

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (needs Python 3 + matplotlib):
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###
set -u

mkdir -p "$OUT_DIR"

python3 -m pangenome_pipeline.step05_cross_genome_annotation_integration.plot_filter_qc \
  --strains "$STRAINS" \
  --original-dir "$ORIGINAL_DIR" \
  --filtered-dir "$FILTERED_DIR" \
  --min-overlap-bp "$MIN_OVERLAP_BP" \
  --out-dir "$OUT_DIR"

echo "QC plots completed. Outputs are in $OUT_DIR"
