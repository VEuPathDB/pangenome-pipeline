#!/bin/bash
#SBATCH --job-name=build_explorer
#SBATCH --partition=low
#SBATCH --output=pangenome_explorer/logs/%x_%j.out
#SBATCH --error=pangenome_explorer/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail

# The interactive pangenome explorer (explorer.html + data/), and the same page packed in one file
# (explorer_standalone.html, opens with a double click, no server). Two uses:
#   preliminary (after the safe fix): the fixed matrix and step 05's GFFs, to look for genes to remove
#       mkdir -p pangenome_explorer/logs && sbatch --export=ALL,SET=preliminary 07_Post_GET_PANGENES_adjusting/scripts/build_explorer_job.sh
#   final (last part of step 07, default): the final pangenome
#       mkdir -p pangenome_explorer/logs && sbatch 07_Post_GET_PANGENES_adjusting/scripts/build_explorer_job.sh
# The output goes to $OUT_DIR/$SET. The log folder must exist BEFORE submitting.

# ===== Edit these =====
SET="${SET-final}"                                           # preliminary | final
TITLE="${TITLE-Toxo Pangenome Explorer}"
FIXED_DIR="${FIXED_DIR-pangenome_fixed}"                     # preliminary: fix_pangenome_safe_job.sh output
STEP05_GFF_DIR="${STEP05_GFF_DIR-basket_annotation/gff_prefer_panel_model_final}"  # preliminary: step 05 final GFFs
FINAL_DIR="${FINAL_DIR-pangenome_final}"                     # final: remove_panel_genes_job.sh output
WINNERS="${WINNERS-basket_annotation/winners.tsv}"           # step 05: PSAURON of every panel gene (optional)
ENTRY_LIST="${ENTRY_LIST-$HOME/hc-storage/interpro_entry.list}"   # InterPro names (optional; IDs only when missing)
OUT_DIR="${OUT_DIR-pangenome_explorer}"
CORE_MISSING="${CORE_MISSING-0}"          # core = missing from at most this many genomes
SOFTCORE_MISSING="${SOFTCORE_MISSING-1}"  # soft-core = missing from at most this many genomes (and not core)
CLOUD_MAX="${CLOUD_MAX-2}"                # cloud = present in at most this many genomes; shell = the rest
# ======================

set +u
module purge 2>/dev/null || true
### ADAPT ### — put python3 on the PATH (standard library only):
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
    [[ -f "$cand/build_explorer_data.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find build_explorer_data.py (run from the project folder)" >&2; exit 1; }

ABSENCES=""
if [[ "$SET" == preliminary ]]; then
    MATRIX="$FIXED_DIR/pangene_matrix_genes.tr.tab"
    GFF_DIR="$STEP05_GFF_DIR"
    [[ -f "$FIXED_DIR/absences_explained.tsv" ]] && ABSENCES="$FIXED_DIR/absences_explained.tsv"
    SUBTITLE="Preliminary pangenome (after the safe fix, before the removal of curated genes)"
elif [[ "$SET" == final ]]; then
    MATRIX="$FINAL_DIR/pangene_matrix_genes.tr.tab"
    GFF_DIR="$FINAL_DIR/gff"
    [[ -f "$FINAL_DIR/absences_explained.tsv" ]] && ABSENCES="$FINAL_DIR/absences_explained.tsv"
    SUBTITLE="Final pangenome (after the removal of curated genes)"
else
    echo "ERROR: SET must be preliminary or final" >&2; exit 1
fi
[[ -f "$MATRIX" ]] || { echo "ERROR: missing $MATRIX" >&2; exit 1; }
[[ -d "$GFF_DIR" ]] || { echo "ERROR: missing folder $GFF_DIR" >&2; exit 1; }
echo "set $SET | matrix $MATRIX | GFFs $GFF_DIR | absences ${ABSENCES:-none}"

O="$OUT_DIR/$SET"
mkdir -p "$O/data"
EXTRA=()
[[ -n "$ABSENCES" ]] && EXTRA+=(--absences "$ABSENCES")
[[ -f "$WINNERS" ]] && EXTRA+=(--winners "$WINNERS")
[[ -f "$ENTRY_LIST" ]] && EXTRA+=(--entry-list "$ENTRY_LIST")
python3 "$SD/build_explorer_data.py" --genes-matrix "$MATRIX" --gff-dir "$GFF_DIR" ${EXTRA[@]+"${EXTRA[@]}"} -o "$O/data" \
  --core-missing "$CORE_MISSING" --softcore-missing "$SOFTCORE_MISSING" --cloud-max "$CLOUD_MAX"
python3 "$SD/build_explorer_locus_index.py" --genes-matrix "$MATRIX" --gff-dir "$GFF_DIR" -o "$O/data"
sed -e "s|__TITLE__|${TITLE}|g" -e "s|__SUBTITLE__|${SUBTITLE}|g" "$SD/explorer.html" > "$O/explorer.html"
cp "$SD/serve_explorer.py" "$O/"
python3 "$SD/make_standalone_explorer.py" "$O" "$O/explorer_standalone.html"

echo "Done. Copy $O/explorer_standalone.html to your computer and open it with a double click."
