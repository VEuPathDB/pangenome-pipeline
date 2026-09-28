#!/bin/bash
#SBATCH --job-name=plot_occupancy_properties
#SBATCH --partition=low
#SBATCH --output=pangenome_plots/logs/%x_%j.out
#SBATCH --error=pangenome_plots/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Run this after get_pangenes has finished. Submit from your project folder. The log folder must
# exist BEFORE submitting (SLURM does not create it):
#     mkdir -p pangenome_plots/logs && sbatch 08_Pangenome_Visualisation/scripts/plot_occupancy_properties_job.sh

# Edit these for your run.
PANGENES_OUT_DIR="for_get_pangenes_7b_pangenes"   # get_pangenes output folder (holds the *2neigh_algMmap_ subfolder)
GFF_DIR="for_get_pangenes_7b"                     # the GFFs given to get_pangenes (link_for_get_pangenes_7b.slurm.sh's OUT_DIR)
GFF_SUFFIX="_for_get_pangenes.gff"                # <genome><suffix> inside GFF_DIR
OUT_DIR="pangenome_plots"
LABEL="${LABEL:?set LABEL (prefix of every output file, e.g. Toxo7b) here or with sbatch --export=ALL,LABEL=...}"

set +u
module purge
source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
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
