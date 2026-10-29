#!/bin/bash
#SBATCH --job-name=inspect_pangenome
#SBATCH --partition=low
#SBATCH --output=pangenome_inspection/logs/%x_%j.out
#SBATCH --error=pangenome_inspection/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

# Step 07, inspection of the pangenome built on step 05's basket annotation. Nothing is changed:
#   1. audit_absences.py              every absence (a pangene with no gene in a genome), and why
#   2. inspect_multigene_pangenes.py  every pangene with more than one gene in a genome, and why
# Submit from your project folder. The log folder must exist BEFORE submitting:
#     mkdir -p pangenome_inspection/logs && sbatch 07_Post_GET_PANGENES_adjusting/scripts/inspect_pangenome_job.sh
# Any setting below can be changed on the command line, e.g.:
#     sbatch --export=ALL,GENOME_DIR=genomes 07_Post_GET_PANGENES_adjusting/scripts/inspect_pangenome_job.sh

# ===== Edit these =====
PANGENES_DIR="${PANGENES_DIR-for_get_pangenes_basket_prefer_panel_model_pangenes}"   # get_pangenes output folder
MATRIX="${MATRIX-}"                                          # or give pangene_matrix_genes.tr.tab directly
GFF_DIR="${GFF_DIR-basket_annotation/gff_prefer_panel_model_final}"   # the GFFs given to get_pangenes (step 05)
GFF_SUFFIX="${GFF_SUFFIX-_basket_annotation.gff}"
GENOMES="${GENOMES-genome_order.txt}"                        # the reference first
GENOME_DIR="${GENOME_DIR-genomes}"                           # <genome>.fna, as given to step 05
PAIRWISE_DIR="${PAIRWISE_DIR-liftoff_output/polished}"       # the Liftoff transfers step 05 read
STEP05_SCRIPTS="${STEP05_SCRIPTS-05_Basket_annotation_integration/scripts}"   # step 05's check_orf.py, decide_models.py
OUT_DIR="${OUT_DIR-pangenome_inspection}"
THREADS="${THREADS-8}"                                       # match --cpus-per-task above
# ======================

set +u
module purge
### ADAPT ### — put blastn / blastp / makeblastdb / python3 on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u
for tool in python3 blastn blastp makeblastdb; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool is not on the PATH (see ### ADAPT ###)" >&2; exit 1; }
done

# Locate the scripts. Under `sbatch` the wrapper is copied to a spool dir, so ${BASH_SOURCE[0]} is useless there -
# fall back to the standard layout relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts" \
    "07_Post_GET_PANGENES_adjusting/scripts"
do
    [[ -f "$cand/inspect_multigene_pangenes.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find inspect_multigene_pangenes.py (run from the project folder)" >&2; exit 1; }
[[ -f "$STEP05_SCRIPTS/check_orf.py" ]] || { echo "ERROR: no check_orf.py in STEP05_SCRIPTS=$STEP05_SCRIPTS" >&2; exit 1; }

if [[ -z "$MATRIX" ]]; then
    for x in "$PANGENES_DIR/pangene_matrix_genes.tr.tab" "$PANGENES_DIR"/*2neigh_algMmap_/pangene_matrix_genes.tr.tab; do
        [[ -f "$x" ]] && { MATRIX="$x"; break; }
    done
    [[ -n "$MATRIX" ]] || { echo "ERROR: no pangene_matrix_genes.tr.tab in $PANGENES_DIR; set MATRIX" >&2; exit 1; }
fi
REF=$(grep -v '^#' "$GENOMES" | awk 'NF{print $1; exit}')
echo "matrix: $MATRIX; reference: $REF"

echo "== 1. absences"
python3 "$SD/audit_absences.py" --genes-matrix "$MATRIX" --gff-dir "$GFF_DIR" --gff-suffix "$GFF_SUFFIX" \
  --genome-dir "$GENOME_DIR" --pairwise-dir "$PAIRWISE_DIR" --threads "$THREADS" --out-dir "$OUT_DIR/absence_audit"
echo "== 2. pangenes with more than one gene in a genome"
python3 "$SD/inspect_multigene_pangenes.py" --genes-matrix "$MATRIX" --gff-dir "$GFF_DIR" --gff-suffix "$GFF_SUFFIX" \
  --genome-dir "$GENOME_DIR" --genomes "$GENOMES" --reference "$REF" --step05-scripts "$STEP05_SCRIPTS" \
  --out-dir "$OUT_DIR/multigene"

echo "Done. $OUT_DIR/absence_audit/summary.tsv and $OUT_DIR/multigene/summary.tsv"
