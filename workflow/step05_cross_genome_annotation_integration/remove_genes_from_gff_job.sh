#!/bin/bash
#SBATCH --job-name=remove_genes_from_gff
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#
# Removes the genes listed in a manual kill list (and everything nested under them) from each
# genome's GFF - see remove_genes_from_gff.py's own header for exactly what it does.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/remove_genes_from_gff_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── Edit these for your run ───
GENOMES="genome_order.txt"                         # comma-separated genome names, or a file with one per line
GFF_DIR="$STEP04_OUT/gff_with_interpro"           # GFFs to remove genes from (never edited in place)
KILL_LIST="$STEP05_OUT/kill_list_manual.tsv"      # written by hand - TSV, no header: genome<TAB>gene_id
OUT_DIR="$STEP05_OUT/gff_curated"                 # where the curated GFFs go (every genome, kill-listed or not)
# ────────────────────────────────

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (standard library only):
activate_env "$ENV_BIOINFO"  # env names / backend: pipeline_config.sh
### END ADAPT ###
set -u

[[ -f "$KILL_LIST" ]] || { echo "ERROR: kill list not found: $KILL_LIST" >&2; exit 1; }

python3 -m pangenome_pipeline.step05_cross_genome_annotation_integration.remove_genes_from_gff \
  --genomes "$GENOMES" --gff-dir "$GFF_DIR" \
  --kill-list "$KILL_LIST" --out-dir "$OUT_DIR"

echo "Done. Curated GFFs are in $OUT_DIR"
