#!/bin/bash
#SBATCH --job-name=explore_multigene_pangenes
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

# Run this after get_pangenes has finished. Submit from your project folder (logs/ must exist
# there - SLURM does not create it):
#     sbatch $PIPELINE_DIR/workflow/step07_post_get_pangenes_adjusting/explore_multigene_pangenes_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
PANGENES_OUT_DIR="$STEP06_OUT/for_get_pangenes_7b_pangenes"   # get_pangenes output folder (holds the *2neigh_algMmap_ subfolder)
TAGGED_DIR="$STEP05_OUT/tagged_annotation"   # step 05's 7b OUT_DIR (tagged GFFs + <genome>_intersect_pairs.tsv)
OUT_DIR="$STEP07_OUT/multigene_pangenes"
LABEL="${LABEL:?set LABEL (prefix of every output file, e.g. Toxo7b) here or with sbatch --export=ALL,LABEL=...}"

set +u
module purge
### ADAPT ### — put python3 on the PATH (standard library only, nothing else needed):
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###
set -u

mkdir -p "$OUT_DIR"

# get_pangenes writes pangene_matrix_genes.tr.tab inside <genome>_..._2neigh_algMmap_/ under its output folder.
shopt -s nullglob
MATRICES=("$PANGENES_OUT_DIR"/*2neigh_algMmap_/pangene_matrix_genes.tr.tab)
shopt -u nullglob
if [[ ${#MATRICES[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one $PANGENES_OUT_DIR/*2neigh_algMmap_/pangene_matrix_genes.tr.tab, found ${#MATRICES[@]}" >&2
    exit 1
fi

python3 -m pangenome_pipeline.step07_post_get_pangenes_adjusting.explore_multigene_pangenes \
  --genes-matrix "${MATRICES[0]}" \
  --tagged-dir "$TAGGED_DIR" \
  -o "$OUT_DIR/${LABEL}_multigene"

echo "Multi-gene pangene exploration completed. Outputs are in $OUT_DIR"
