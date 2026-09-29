#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=rename_headers
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Trim FASTA headers, in place, in genomes/ (see docs/step01_genome_selection_qc.md).
# Submit from the project folder:
#     sbatch $PIPELINE_DIR/workflow/step01_genome_selection_qc/rename_headers_slurm.sh
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#  1. The  #SBATCH  lines above (partition / resources).
#  2. The  module load  line below (### ADAPT ###) — it only needs a  python3 .
#     Use your cluster's python module, or  source activate <env>, or delete the
#     line if python3 is already on the PATH.
# ═══════════════════════════════════════════════════════════════════════════════

set -eo pipefail

### ADAPT ### — put python3 on the PATH:
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###

set -u

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# --- edit these if your folders differ -------------------------------------
GENOME_DIR="$PROJECT_DIR/genomes"
GFF_DIR="$PROJECT_DIR/gff_originals"     # set to "" to skip the GFF ID cross-check
# --------------------------------------------------------------------------


shopt -s nullglob
inputs=( "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fa )
shopt -u nullglob
[[ ${#inputs[@]} -ge 1 ]] || { echo "ERROR: no FASTA files in $GENOME_DIR" >&2; exit 1; }

python3 -m pangenome_pipeline.step01_genome_selection_qc.rename_headers "${inputs[@]}" --in-place

if [[ -n "$GFF_DIR" ]]; then
    shopt -s nullglob
    fastas=( "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fa )
    shopt -u nullglob
    python3 -m pangenome_pipeline.step01_genome_selection_qc.rename_headers "${fastas[@]}" --check-gff "$GFF_DIR"
fi

echo "Done. $GENOME_DIR/ headers are trimmed."
