#!/bin/bash
#SBATCH --job-name=analyze_input_gffs
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#
# Profiles every genome's own input GFF before Liftoff runs - see analyze_input_gffs.py's own
# header for exactly what it checks and what it writes.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/analyze_input_gffs_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── Edit these for your run ───
GENOMES="genome_order.txt"   # comma-separated genome names, or a file with one per line
GFF_DIR="$STEP04_OUT/gff_with_interpro"   # step 04 output (domain evidence written in) - run this check once, after InterProScan
OUT_DIR="$STEP05_OUT/annotation_check"   # where to write the tables and plots
# ────────────────────────────────

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (needs Python 3 + matplotlib):
activate_env "$ENV_SYNTENY"  # env names / micromamba: pipeline_config.sh
### END ADAPT ###
set -u

python3 -m pangenome_pipeline.step05_cross_genome_annotation_integration.analyze_input_gffs \
  --genomes "$GENOMES" --gff-dir "$GFF_DIR" --out-dir "$OUT_DIR"

echo "Annotation check completed. Outputs are in $OUT_DIR"
