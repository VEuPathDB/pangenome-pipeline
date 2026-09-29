#!/bin/bash
#SBATCH --job-name=prepare_for_liftoff
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#
# Prepares each genome's GFF for Liftoff: unifies gene feature types (so Liftoff transfers the
# protein-coding genes and not the non-coding RNA genes) and removes UTRs. See prepare_for_liftoff.py.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/prepare_for_liftoff_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── Edit these for your run ───
GENOMES="genome_order.txt"                         # comma-separated genome names, or a file with one per line
GFF_DIR="$STEP05_OUT/gff_curated"                 # GFFs to prepare (never edited in place) - use $STEP04_OUT/gff_with_interpro if no genes were kill-listed
OUT_DIR="$STEP05_OUT/gff_for_liftoff"             # where the prepared GFFs (+ liftoff_feature_types.txt) go
# ────────────────────────────────

set +u
module purge
### ADAPT ### — put python3 on the PATH (standard library only, nothing else needed):
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###
set -u

python3 -m pangenome_pipeline.step05_cross_genome_annotation_integration.prepare_for_liftoff \
  --genomes "$GENOMES" --gff-dir "$GFF_DIR" \
  --out-dir "$OUT_DIR"

echo "Done. Prepared GFFs are in $OUT_DIR"
