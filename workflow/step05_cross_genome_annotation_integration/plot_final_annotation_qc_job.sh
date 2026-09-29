#!/bin/bash
#SBATCH --job-name=plot_final_annotation_qc
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

# Submit from your project folder (logs/ must exist there - SLURM does not create it):
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/plot_final_annotation_qc_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
STRAINS="genome_order.txt"                   # comma-separated names or a file with one per line, first = the reference (same as ORDER)
FINAL_DIR="$STEP05_OUT/expanded_annotation"   # build_final_annotation_basic.sh's FINAL_DIR
GENOME_DIR="$STEP05_OUT/gff_for_liftoff"   # each genome's own annotation, to compare the gene counts with
GENOME_SUFFIX=".gff3"                        # filename suffix for each genome's own GFF
OUT_DIR="$STEP05_OUT/final_annotation_qc_plots"   # where to write the plots and summary.tsv

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (needs Python 3 + matplotlib):
activate_env "$ENV_SYNTENY"  # env names / micromamba: pipeline_config.sh
### END ADAPT ###
set -u

python3 -m pangenome_pipeline.step05_cross_genome_annotation_integration.plot_final_annotation_qc \
  --strains "$STRAINS" --final-dir "$FINAL_DIR" \
  --genome-dir "$GENOME_DIR" --genome-suffix "$GENOME_SUFFIX" \
  --out-dir "$OUT_DIR"

echo "QC plots completed. Outputs are in $OUT_DIR"
