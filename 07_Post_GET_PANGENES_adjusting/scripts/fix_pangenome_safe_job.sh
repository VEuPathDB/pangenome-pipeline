#!/bin/bash
#SBATCH --job-name=fix_pangenome_safe
#SBATCH --partition=low
#SBATCH --output=pangenome_fixed/logs/%x_%j.out
#SBATCH --error=pangenome_fixed/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Step 07, after inspect_pangenome_job.sh: apply only the safe fixes to the pangene matrix (genes that get_pangenes
# split into two pangenes are moved back together) and explain every absence left. The GFFs are not changed.
# Then summarize_pangenes.py lists, on the fixed pangenome, the low-confidence pangenes (invalid ORF in more than 50% of
# the strains, with the defect, whether the strains share it, PSAURON and domains) and the genes in no pangene, to
# help decide what to remove.
# Submit from your project folder. The log folder must exist BEFORE submitting:
#     mkdir -p pangenome_fixed/logs && sbatch 07_Post_GET_PANGENES_adjusting/scripts/fix_pangenome_safe_job.sh

# ===== Edit these =====
PANGENES_DIR="${PANGENES_DIR-for_get_pangenes_basket_prefer_panel_model_pangenes}"   # get_pangenes output folder
MATRIX="${MATRIX-}"                                          # or give pangene_matrix_genes.tr.tab directly
ABSENCES="${ABSENCES-pangenome_inspection/absence_audit/absences.tsv}"   # from inspect_pangenome_job.sh
GFF_DIR="${GFF_DIR-basket_annotation/gff_prefer_panel_model_final}"     # the GFFs given to get_pangenes (step 05)
ABSENT="${ABSENT-basket_annotation/decide/absent.tsv}"                  # step 05
STACKED="${STACKED-basket_annotation/gff_prefer_panel_model/stacked_dropped.tsv}"   # step 05
WINNERS="${WINNERS-basket_annotation/winners.tsv}"           # step 05: PSAURON of every panel gene
GENOME_DIR="${GENOME_DIR-genomes}"                           # <genome>.fna, as given to step 05 (where invalid ORFs break)
ENTRY_LIST="${ENTRY_LIST-$HOME/hc-storage/interpro_entry.list}"   # InterPro names (optional)
OUT_DIR="${OUT_DIR-pangenome_fixed}"
# ======================

set +u
module purge
### ADAPT ### — put python3 on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

# Locate the script. Under `sbatch` the wrapper is copied to a spool dir, so ${BASH_SOURCE[0]} is useless there -
# fall back to the standard layout relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts" \
    "07_Post_GET_PANGENES_adjusting/scripts"
do
    [[ -f "$cand/fix_pangenome_safe.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find fix_pangenome_safe.py (run from the project folder)" >&2; exit 1; }
for f in "$ABSENCES" "$ABSENT" "$STACKED"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

if [[ -z "$MATRIX" ]]; then
    for x in "$PANGENES_DIR/pangene_matrix_genes.tr.tab" "$PANGENES_DIR"/*2neigh_algMmap_/pangene_matrix_genes.tr.tab; do
        [[ -f "$x" ]] && { MATRIX="$x"; break; }
    done
    [[ -n "$MATRIX" ]] || { echo "ERROR: no pangene_matrix_genes.tr.tab in $PANGENES_DIR; set MATRIX" >&2; exit 1; }
fi
echo "matrix: $MATRIX"

python3 "$SD/fix_pangenome_safe.py" --genes-matrix "$MATRIX" --absences "$ABSENCES" --gff-dir "$GFF_DIR" \
  --absent "$ABSENT" --stacked-dropped "$STACKED" --out-dir "$OUT_DIR"

echo "== lists on the fixed pangenome: low-confidence pangenes, genes in no pangene"
EXTRA=()
[[ -f "$ENTRY_LIST" ]] && EXTRA+=(--entry-list "$ENTRY_LIST")
PYTHONPATH="$SD" python3 "$SD/summarize_pangenes.py" --genes-matrix "$OUT_DIR/pangene_matrix_genes.tr.tab" \
  --gff-dir "$GFF_DIR" --winners "$WINNERS" --genome-dir "$GENOME_DIR" ${EXTRA[@]+"${EXTRA[@]}"} \
  --out-dir "$OUT_DIR/pangene_lists"

echo "Done. $OUT_DIR/summary.tsv, moves.tsv, absences_explained.tsv, pangene_lists/"
