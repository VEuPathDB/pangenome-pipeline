#!/bin/bash
#SBATCH --job-name=plot_basket_annotation
#SBATCH --partition=low
#SBATCH --output=basket_annotation/logs/%x_%j.out
#SBATCH --error=basket_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Plots of the basket annotation again, from the final GFFs that basket_annotation_job.sh already wrote (nothing else
# is run), e.g. for one mode only. They go to $OUT_DIR/plots/ named <PLOT_LABEL>_<modes>_*.png (e.g.
# Toxo_prefer_panel_model_orf_quality.png), so the job's own plots are not overwritten.
# Submit from your project folder:
#     sbatch --export=ALL,MODES=prefer_panel_model 05_Basket_annotation_integration/scripts/plot_basket_annotation_job.sh
#     sbatch 05_Basket_annotation_integration/scripts/plot_basket_annotation_job.sh     # the MODES set below

# ===== Edit these (the same values as in basket_annotation_job.sh) =====
GENOMES="genome_order.txt"
OUT_DIR="basket_annotation"
NATIVE_DIR="gff_for_liftoff"                        # used only if $OUT_DIR/input/native (written by the job) is missing
PLOT_LABEL="panel"
MODES="${MODES:-prefer_panel_model}"   # or set with --export=ALL,MODES=... (e.g. "keep_similar prefer_panel_model")
# =======================================================================

set +u
module purge
### ADAPT ### — put python3 with matplotlib on the PATH:
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
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Basket_annotation_integration/scripts" \
    "05_Basket_annotation_integration/scripts"
do
    [[ -f "$cand/plot_basket_annotation.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find plot_basket_annotation.py (run from the project folder)" >&2; exit 1; }
[[ -d "$OUT_DIR/input/native" ]] && NATIVE_DIR="$OUT_DIR/input/native"
for m in $MODES; do
    [[ -d "$OUT_DIR/gff_${m}_final" ]] || { echo "ERROR: $OUT_DIR/gff_${m}_final not found: run basket_annotation_job.sh for mode $m" >&2; exit 1; }
done

LABEL="${PLOT_LABEL}_$(echo $MODES | tr ' ' '_')"   # the modes in the file names: the job's own plots are not overwritten
echo "plots for: $MODES -> $OUT_DIR/plots/${LABEL}_*.png"
python3 "$SD/plot_basket_annotation.py" --result-dir "$OUT_DIR" --genomes "$GENOMES" --native-dir "$NATIVE_DIR" \
  --modes "$MODES" --label "$LABEL" --out-dir "$OUT_DIR/plots"
