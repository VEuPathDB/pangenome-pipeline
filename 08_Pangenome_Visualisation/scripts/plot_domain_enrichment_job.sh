#!/bin/bash
#SBATCH --job-name=plot_domain_enrichment
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
#     mkdir -p pangenome_plots/logs && sbatch 08_Pangenome_Visualisation/scripts/plot_domain_enrichment_job.sh

# Edit these for your run.
PANGENES_OUT_DIR="for_get_pangenes_7b_pangenes"                             # get_pangenes output folder (holds the *2neigh_algMmap_ subfolder)
GFF_DIR="for_get_pangenes_7b"                                             # the GFFs given to get_pangenes (link_for_get_pangenes_7b.slurm.sh's OUT_DIR)
ENTRY_LIST=""                                                             # InterPro IPR id -> name; empty = INTERPRO_ENTRY_LIST from pipeline_config.sh
OUT_DIR="pangenome_plots"
LABEL="${LABEL:?set LABEL (prefix of every output file, e.g. Toxo7b) here or with sbatch --export=ALL,LABEL=...}"
CATEGORIES="Core,Cloud"                                                   # comma-separated: Core, Soft-core, Shell, Cloud
ALPHA="0.05"                                                              # adjusted p-value cutoff
MIN_TOTAL="5"                                                             # entries carried by fewer clusters in the whole universe are not tested
MIN_COUNT="3"                                                             # min clusters carrying the entry in the category, to be drawn
TOP="15"                                                                  # max entries shown per category

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

find_script PY 08_Pangenome_Visualisation/scripts/plot_domain_enrichment.py
mkdir -p "$OUT_DIR"

# get_pangenes writes pangene_matrix_genes.tr.tab inside <genome>_..._2neigh_algMmap_/ under its output folder.
shopt -s nullglob
MATRICES=("$PANGENES_OUT_DIR"/*2neigh_algMmap_/pangene_matrix_genes.tr.tab)
shopt -u nullglob
if [[ ${#MATRICES[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one $PANGENES_OUT_DIR/*2neigh_algMmap_/pangene_matrix_genes.tr.tab, found ${#MATRICES[@]}" >&2
    exit 1
fi
GENES_MATRIX="${MATRICES[0]}"

# Compute nodes may not reach the internet: download the list once on the login node (the curl below),
# otherwise the plot falls back to IPR ids.
[[ -n "$ENTRY_LIST" ]] || ENTRY_LIST="${INTERPRO_ENTRY_LIST:-}"
if [[ -z "$ENTRY_LIST" ]]; then
    echo "WARNING: no InterPro entry list set (INTERPRO_ENTRY_LIST in pipeline_config.sh); terms will be labelled by IPR id" >&2
elif [[ ! -f "$ENTRY_LIST" ]]; then
    curl -sfL -o "$ENTRY_LIST" "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list" \
        || { echo "WARNING: could not download the InterPro entry list; terms will be labelled by IPR id" >&2; rm -f "$ENTRY_LIST"; }
fi
ENTRY_ARGS=()
[[ -n "$ENTRY_LIST" && -f "$ENTRY_LIST" ]] && ENTRY_ARGS+=(--entry-list "$ENTRY_LIST")

python3 "$PY" \
  --genes-matrix "$GENES_MATRIX" --gff-dir "$GFF_DIR" \
  -o "$OUT_DIR/${LABEL}_domain_enrichment" \
  --categories "$CATEGORIES" --alpha "$ALPHA" --min-total "$MIN_TOTAL" --min-count "$MIN_COUNT" --top "$TOP" \
  ${ENTRY_ARGS[@]+"${ENTRY_ARGS[@]}"}

echo "Domain enrichment plot completed. Outputs are in $OUT_DIR"
