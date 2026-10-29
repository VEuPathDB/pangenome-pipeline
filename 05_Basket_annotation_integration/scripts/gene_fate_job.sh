#!/bin/bash
#SBATCH --job-name=gene_fate
#SBATCH --partition=low
#SBATCH --output=basket_annotation/logs/%x_%j.out
#SBATCH --error=basket_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Points 9-10 of the basket annotation alone, on results already written by basket_annotation_job.sh (nothing else is
# run): the fate of every gene (gene_fate.py) and the figures made from it (plot_gene_fate.py), per mode in MODES.
# Submit from your project folder; the log folder must exist BEFORE submitting:
#     sbatch 05_Basket_annotation_integration/scripts/gene_fate_job.sh
# With another results folder, the logs go to basket_annotation/logs all the same (or set --output/--error):
#     sbatch --export=ALL,OUT_DIR=basket_annotation_keep_all_genes 05_Basket_annotation_integration/scripts/gene_fate_job.sh

# ===== Edit these for your species and run =====  (each can also be set on the command line: sbatch --export=ALL,NAME=value ...)
GENOMES="${GENOMES-genome_order.txt}"                           # one genome per line, the REFERENCE FIRST
OUT_DIR="${OUT_DIR-basket_annotation}"                          # the results of basket_annotation_job.sh
MODES="${MODES-prefer_panel_model}"                             # the modes to do (their gff_<mode> and gff_<mode>_final must exist)
PLOT_LABEL="${PLOT_LABEL-panel}"                                # prefix of the plot files
MIN_AA="${MIN_AA-100}"                                          # the same values as in basket_annotation_job.sh
PSAURON_MIN="${PSAURON_MIN-0.7}"
PSAURON_MIN_IN_REF="${PSAURON_MIN_IN_REF-0.9}"
# genes removed before step 05, for the gene balance (optional; several entries separated by "|", each FILE or
# FILE:REASON; see plot_gene_fate.py):
REMOVED_GENES="${REMOVED_GENES-}"                               # tables of genes (genome, gene), or {genome} GFFs, e.g. step 02's TE genes
REMOVED_CONTIGS="${REMOVED_CONTIGS-}"                           # tables of contigs (genome, contig) removed in step 01
ORIGINAL_GFF_DIR="${ORIGINAL_GFF_DIR-}"                         # the GFFs before the contigs were removed (the genes on REMOVED_CONTIGS)
# ===============================================

set +u
module purge
### ADAPT ### — put python3 (with matplotlib and numpy) on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

# Locate the scripts. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there - fall back to the standard layout
# relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Basket_annotation_integration/scripts" \
    "05_Basket_annotation_integration/scripts"
do
    [[ -f "$cand/gene_fate.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find gene_fate.py (run from the project folder)" >&2; exit 1; }
REF=$(grep -v '^#' "$GENOMES" | awk 'NF{print $1; exit}')
[[ -d "$OUT_DIR/input/native" ]] || { echo "ERROR: no $OUT_DIR/input/native (the native GFFs basket_annotation_job.sh read)" >&2; exit 1; }

REMOVED_ARGS=()
IFS='|' read -r -a _rg <<< "$REMOVED_GENES"
for x in ${_rg[@]+"${_rg[@]}"}; do [[ -n "$x" ]] && REMOVED_ARGS+=(--removed-genes "$x"); done
IFS='|' read -r -a _rc <<< "$REMOVED_CONTIGS"
for x in ${_rc[@]+"${_rc[@]}"}; do [[ -n "$x" ]] && REMOVED_ARGS+=(--removed-contigs "$x"); done
[[ -n "$ORIGINAL_GFF_DIR" ]] && REMOVED_ARGS+=(--original-gff-dir "$ORIGINAL_GFF_DIR")

mkdir -p "$OUT_DIR/plots"
for MODE in $MODES; do
    echo "== gene fate, mode $MODE"
    python3 "$SD/gene_fate.py" --native-dir "$OUT_DIR/input/native" --genomes "$GENOMES" --result-dir "$OUT_DIR" \
      --build-dir "$OUT_DIR/gff_$MODE" --final-dir "$OUT_DIR/gff_${MODE}_final" \
      --out "$OUT_DIR/gene_fate_$MODE.tsv"
    echo "== gene balance and eliminated genes, mode $MODE"
    python3 "$SD/plot_gene_fate.py" --gene-fate "$OUT_DIR/gene_fate_$MODE.tsv" \
      --final-dir "$OUT_DIR/gff_${MODE}_final" --result-dir "$OUT_DIR" --reference "$REF" \
      --label "${PLOT_LABEL}_${MODE}" --min-aa "$MIN_AA" --psauron-min "$PSAURON_MIN" \
      --psauron-min-in-ref "$PSAURON_MIN_IN_REF" ${REMOVED_ARGS[@]+"${REMOVED_ARGS[@]}"} --out-dir "$OUT_DIR/plots"
    python3 "$SD/plot_eliminated_psauron.py" --eliminated "$OUT_DIR/plots/${PLOT_LABEL}_${MODE}_eliminated.tsv" \
      --label "${PLOT_LABEL}_${MODE}" --psauron-min "$PSAURON_MIN" --psauron-min-in-ref "$PSAURON_MIN_IN_REF" \
      --out-dir "$OUT_DIR/plots"
done
if grep -q $'\tNO RECORD\t' "$OUT_DIR"/gene_fate_*.tsv; then
    echo "ERROR: native genes without a final model and without a record (fate NO RECORD in $OUT_DIR/gene_fate_<mode>.tsv)" >&2
    exit 1
fi
echo "Done: $OUT_DIR/gene_fate_<mode>.tsv and $OUT_DIR/plots/${PLOT_LABEL}_<mode>_gene_balance.png, _eliminated_reasons.png"
