#!/bin/bash
#SBATCH --job-name=plot_all
#SBATCH --partition=low
#SBATCH --output=pangenome_plots/logs/%x_%j.out
#SBATCH --error=pangenome_plots/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Step 07: the pangenome plots (plot_all.sh) and the InterPro enrichment per occupancy category
# (plot_domain_enrichment.py). Two uses:
#   preliminary (after the safe fix): the fixed matrix and step 05's GFFs, to look for genes to remove; the POCS
#       heatmap comes from get_pangenes' POCS.matrix.tab (the genomes do not change)
#       mkdir -p pangenome_plots/logs && sbatch --export=ALL,SET=preliminary 07_Post_GET_PANGENES_adjusting/scripts/plot_all_job.sh
#   final (last part of step 07, default): the final pangenome (remove_panel_genes_job.sh output)
#       mkdir -p pangenome_plots/logs && sbatch 07_Post_GET_PANGENES_adjusting/scripts/plot_all_job.sh
# The figures go to $OUT_DIR/$SET, every file starts with <LABEL>_<SET>. The log folder must exist BEFORE submitting.

# ===== Edit these =====
SET="${SET-final}"                                           # preliminary | final
LABEL="${LABEL-Toxo}"                                        # species label, e.g. Toxo or Afum
FIXED_DIR="${FIXED_DIR-pangenome_fixed}"                     # preliminary: fix_pangenome_safe_job.sh output
STEP05_GFF_DIR="${STEP05_GFF_DIR-basket_annotation/gff_prefer_panel_model_final}"  # preliminary: step 05 final GFFs
PANGENES_DIR="${PANGENES_DIR-for_get_pangenes_basket_prefer_panel_model_pangenes}"   # preliminary: get_pangenes' POCS
FINAL_DIR="${FINAL_DIR-pangenome_final}"                     # final: remove_panel_genes_job.sh output
ENTRY_LIST="${ENTRY_LIST-$HOME/hc-storage/interpro_entry.list}"   # InterPro id -> name; downloaded here if missing
OUT_DIR="${OUT_DIR-pangenome_plots}"
CORE_MISSING="${CORE_MISSING-0}"          # core = missing from at most this many genomes
SOFTCORE_MISSING="${SOFTCORE_MISSING-1}"  # soft-core = missing from at most this many genomes (and not core)
CLOUD_MAX="${CLOUD_MAX-2}"                # cloud = present in at most this many genomes; shell = the rest
# ======================

set +u
module purge 2>/dev/null || true
### ADAPT ### — put python3 with numpy, scipy and matplotlib on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u
for mod in numpy scipy matplotlib; do
    python3 -c "import $mod" 2>/dev/null || { echo "ERROR: python module '$mod' missing (see ### ADAPT ###)" >&2; exit 1; }
done

# Locate the scripts. Under `sbatch` the wrapper is copied to a spool dir, so ${BASH_SOURCE[0]} is useless there -
# fall back to the standard layout relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts" \
    "07_Post_GET_PANGENES_adjusting/scripts"
do
    [[ -f "$cand/plot_all.sh" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find plot_all.sh (run from the project folder)" >&2; exit 1; }

POCS=""
if [[ "$SET" == preliminary ]]; then
    MATRIX="$FIXED_DIR/pangene_matrix_genes.tr.tab"
    GENES_MATRIX="$MATRIX"
    GFF_DIR="$STEP05_GFF_DIR"
    # get_pangenes writes POCS.matrix.tab in a subfolder (<genome>_..._2neigh_algMmap_)
    for x in "$PANGENES_DIR/POCS.matrix.tab" "$PANGENES_DIR"/*2neigh_algMmap_/POCS.matrix.tab; do
        [[ -f "$x" ]] && { POCS="$x"; break; }
    done
elif [[ "$SET" == final ]]; then
    MATRIX="$FINAL_DIR/pangene_matrix_genes.tr.tab"     # the POCS heatmap is made with the preliminary set only
    GENES_MATRIX="$MATRIX"
    GFF_DIR="$FINAL_DIR/gff"
else
    echo "ERROR: SET must be preliminary or final" >&2; exit 1
fi
for f in "$MATRIX" "$GENES_MATRIX"; do [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }; done
echo "set $SET | matrix $MATRIX | GFFs $GFF_DIR"

# Compute nodes may not reach the internet: the list is downloaded once if missing; without it the terms are IPR ids
if [[ ! -f "$ENTRY_LIST" ]]; then
    curl -sfL -o "$ENTRY_LIST" "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list" \
        || { echo "WARNING: no InterPro entry list; terms will be labelled by IPR id" >&2; rm -f "$ENTRY_LIST"; }
fi
ENTRY_ARGS=()
[[ -f "$ENTRY_LIST" ]] && ENTRY_ARGS+=(--entry-list "$ENTRY_LIST")

O="$OUT_DIR/$SET"
export CORE_MISSING SOFTCORE_MISSING CLOUD_MAX      # read by plot_all.sh
bash "$SD/plot_all.sh" -i "$MATRIX" -o "$O" -l "${LABEL}_${SET}"
python3 "$SD/plot_domain_enrichment.py" --genes-matrix "$GENES_MATRIX" --gff-dir "$GFF_DIR" \
    -o "$O/${LABEL}_${SET}_domain_enrichment" ${ENTRY_ARGS[@]+"${ENTRY_ARGS[@]}"} \
    --core-missing "$CORE_MISSING" --softcore-missing "$SOFTCORE_MISSING" --cloud-max "$CLOUD_MAX"
[[ -n "$POCS" ]] && python3 "$SD/plot_pocs.py" -i "$POCS" -o "$O/${LABEL}_${SET}_POCS_heatmap"
find "$O" -name "*.svg" -delete          # PNG only

echo "Done. Figures in $O"
