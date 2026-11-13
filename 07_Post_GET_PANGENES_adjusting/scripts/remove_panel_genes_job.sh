#!/bin/bash
#SBATCH --job-name=remove_panel_genes
#SBATCH --partition=low
#SBATCH --output=pangenome_final/logs/%x_%j.out
#SBATCH --error=pangenome_final/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Step 07, after fix_pangenome_safe_job.sh: remove from the final GFFs and from the matrix the whole panel genes
# (e.g. transposable elements) listed in REMOVE_LIST (columns panel_gene, reason), and the single genes listed in
# REMOVE_GENE_LIST (columns gene, reason; the final ID or the genome's own gene ID). Empty lists (header only)
# copy everything unchanged, so pangenome_final/ always holds the result of step 07. Then final_tables.py writes
# gene_status.tsv (every gene: kept, moved, not in a pangene or removed, and why) and singletons.tsv (why each
# single-strain pangene is alone). Add genes to the list and rerun as many times as needed: it always starts from FIXED_DIR.
# Submit from your project folder. The log folder must exist BEFORE submitting:
#     mkdir -p pangenome_final/logs && sbatch 07_Post_GET_PANGENES_adjusting/scripts/remove_panel_genes_job.sh

# ===== Edit these =====
REMOVE_LIST="${REMOVE_LIST-07_Post_GET_PANGENES_adjusting/genes_to_remove.tsv}"                  # whole panel genes
REMOVE_GENE_LIST="${REMOVE_GENE_LIST-07_Post_GET_PANGENES_adjusting/strain_genes_to_remove.tsv}"   # single genes
EXCLUDE_LIST="${EXCLUDE_LIST-07_Post_GET_PANGENES_adjusting/genes_to_exclude_from_pangenome.tsv}"   # panel genes left out of the pangenome only (used if present)
ACTIONS="${ACTIONS-07_Post_GET_PANGENES_adjusting/curation_actions.tsv}"   # the reference check page's table (used if present)
GFF_DIR="${GFF_DIR-basket_annotation/gff_prefer_panel_model_final}"   # step 05 final GFFs
FIXED_DIR="${FIXED_DIR-pangenome_fixed}"                               # fix_pangenome_safe_job.sh output
PANGENES_DIR="${PANGENES_DIR-for_get_pangenes_basket_prefer_panel_model_pangenes}"   # get_pangenes output folder
STEP05_DIR="${STEP05_DIR-basket_annotation}"                 # step 05 (winners.tsv, not_eligible.tsv)
PANGENE_PREFIX="${PANGENE_PREFIX-}"                       # final pangene IDs <prefix>_pan_00001... (Tg, Af); "" = keep get_pangenes' names
OUT_DIR="${OUT_DIR-pangenome_final}"
# ======================

set +u
module purge 2>/dev/null || true
### ADAPT ### — put python3 on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

# Locate the scripts. Under `sbatch` the wrapper is copied to a spool dir, so ${BASH_SOURCE[0]} is useless there -
# fall back to the standard layout relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts" \
    "07_Post_GET_PANGENES_adjusting/scripts"
do
    [[ -f "$cand/remove_panel_genes.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find remove_panel_genes.py (run from the project folder)" >&2; exit 1; }
for f in "$REMOVE_LIST" "$REMOVE_GENE_LIST" "$FIXED_DIR/pangene_matrix_genes.tr.tab" "$FIXED_DIR/absences_explained.tsv"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

AO=()
if [[ -n "$ACTIONS" && -f "$ACTIONS" ]]; then AO=(--actions "$ACTIONS"); echo "curation actions: $ACTIONS"
else echo "curation actions: none ($ACTIONS not found)"; fi
if [[ -n "$EXCLUDE_LIST" && -f "$EXCLUDE_LIST" ]]; then AO+=(--exclude "$EXCLUDE_LIST"); echo "left out of the pangenome: $EXCLUDE_LIST"
else echo "left out of the pangenome: none ($EXCLUDE_LIST not found)"; fi
PYTHONPATH="$SD" python3 "$SD/remove_panel_genes.py" --remove "$REMOVE_LIST" --remove-genes "$REMOVE_GENE_LIST" ${AO[@]+"${AO[@]}"} \
  --gff-dir "$GFF_DIR" \
  --genes-matrix "$FIXED_DIR/pangene_matrix_genes.tr.tab" --absences "$FIXED_DIR/absences_explained.tsv" \
  --pangene-prefix "$PANGENE_PREFIX" --out-dir "$OUT_DIR"

MATRIX=""
for x in "$PANGENES_DIR/pangene_matrix_genes.tr.tab" "$PANGENES_DIR"/*2neigh_algMmap_/pangene_matrix_genes.tr.tab; do
    [[ -f "$x" ]] && { MATRIX="$x"; break; }
done
[[ -n "$MATRIX" ]] || { echo "ERROR: no pangene_matrix_genes.tr.tab in $PANGENES_DIR" >&2; exit 1; }
PYTHONPATH="$SD" python3 "$SD/final_tables.py" --step05-gff-dir "$GFF_DIR" --step05-dir "$STEP05_DIR" \
  --pangenes-matrix "$MATRIX" --moves "$FIXED_DIR/moves.tsv" --final-dir "$OUT_DIR"

echo "Done. $OUT_DIR/gff, $OUT_DIR/pangene_matrix_genes.tr.tab, pangene_ids.tsv, removed_genes.tsv, excluded_from_pangenome.tsv, genes_not_found.tsv, gene_status.tsv, singletons.tsv"
