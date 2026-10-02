#!/bin/bash
#SBATCH --job-name=TE_summary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Summarize TE content across all genomes from the RepeatMasker output folders
# (RM_<genome_id>/, produced by run_repeatmasker.slurm.sh) and write the summary
# TSVs + composition plots. Light job — a couple of cores and a few GB are plenty.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/workflow/step02_te_annotation_clean_up/summarize_te_content.slurm.sh
#
set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
INDIR="$STEP02_OUT/repeatmasker"   # folder holding the RM_<genome_id>/ dirs (each with one .out + one .tbl)
OUTDIR="$STEP02_OUT/te_content_summary"   # where the TSVs, filtered .out files and plots/ go
TOP_SUPERFAMILIES=8              # superfamilies shown individually before folding into "Other"
TOP_FAMILIES=15                  # families shown individually before folding into "Other"
# ────────────────────────────────────────────────────────────────────────────

### ADAPT ### — put a python3 with pandas + matplotlib on the PATH
#              (e.g. the 'synteny' env from step 01 already has matplotlib;
#               `micromamba install -c conda-forge pandas matplotlib` otherwise):
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
python3 -c "import pandas, matplotlib" 2>/dev/null \
    || { echo "ERROR: python3 is missing pandas and/or matplotlib" >&2; exit 1; }
compgen -G "$INDIR"/RM_*/ >/dev/null || compgen -G "$INDIR"/R_*/ >/dev/null \
    || { echo "ERROR: no RM_<genome_id>/ (or R_<genome_id>/) folders in $INDIR (run run_repeatmasker.slurm.sh first)" >&2; exit 1; }

python3 -m pangenome_pipeline.step02_te_annotation_clean_up.summarize_te_content \
    --indir "$INDIR" \
    --outdir "$OUTDIR" \
    --top-superfamilies "$TOP_SUPERFAMILIES" \
    --top-families "$TOP_FAMILIES"

echo "Done — summary and plots in $OUTDIR/"
