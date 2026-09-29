#!/bin/bash
#SBATCH --job-name=plot_occupancy_properties
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

# Run this after get_pangenes has finished. Submit from your project folder (logs/ must exist
# there - SLURM does not create it):
#     sbatch $PIPELINE_DIR/workflow/step08_pangenome_visualisation/plot_occupancy_properties_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
PANGENES_OUT_DIR="$STEP06_OUT/for_get_pangenes_7b_pangenes"   # get_pangenes output folder (holds the *2neigh_algMmap_ subfolder)
GFF_DIR="$STEP06_OUT/for_get_pangenes_7b"                     # the GFFs given to get_pangenes (link_for_get_pangenes_7b.slurm.sh's OUT_DIR)
GFF_SUFFIX="_for_get_pangenes.gff"                # <genome><suffix> inside GFF_DIR
OUT_DIR="$STEP08_OUT/pangenome_plots"
LABEL="${LABEL:?set LABEL (prefix of every output file, e.g. Toxo7b) here or with sbatch --export=ALL,LABEL=...}"

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (needs Python 3 + numpy + scipy + matplotlib):
activate_env "$ENV_SYNTENY"  # env names / micromamba: pipeline_config.sh
### END ADAPT ###
set -u

for mod in numpy scipy matplotlib; do
    python3 -c "import $mod" 2>/dev/null \
        || { echo "ERROR: the activated env has no python module '$mod' (needs numpy, scipy and matplotlib)" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"

# get_pangenes writes pangene_matrix_genes.tr.tab inside <genome>_..._2neigh_algMmap_/ under its output folder.
shopt -s nullglob
MATRICES=("$PANGENES_OUT_DIR"/*2neigh_algMmap_/pangene_matrix_genes.tr.tab)
shopt -u nullglob
if [[ ${#MATRICES[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one $PANGENES_OUT_DIR/*2neigh_algMmap_/pangene_matrix_genes.tr.tab, found ${#MATRICES[@]}" >&2
    exit 1
fi

# Two versions of the protein-length panel: linear axis and log axis.
python3 -m pangenome_pipeline.step08_pangenome_visualisation.plot_occupancy_properties \
  --genes-matrix "${MATRICES[0]}" --gff-dir "$GFF_DIR" --gff-suffix "$GFF_SUFFIX" \
  -o "$OUT_DIR/${LABEL}_occupancy_properties"
python3 -m pangenome_pipeline.step08_pangenome_visualisation.plot_occupancy_properties \
  --genes-matrix "${MATRICES[0]}" --gff-dir "$GFF_DIR" --gff-suffix "$GFF_SUFFIX" --yscale log \
  -o "$OUT_DIR/${LABEL}_occupancy_properties_log"

echo "Occupancy properties plot completed. Outputs are in $OUT_DIR"
