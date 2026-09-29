#!/bin/bash
#SBATCH --job-name=fix_rm_contigs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Afum-only fix. RepeatMasker was run on the assemblies BEFORE the step-01 contig
# drop, so each RM_<genome_id>/ still has hits on contigs that were removed after.
# This rewrites the .out / .gff / .tbl to match the final dataset, into a new
# folder (originals untouched). Pure Python, stdlib only — a few seconds.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/species_specific/Afum/step02_te_annotation_clean_up/fix_repeatmasker_after_contig_drop.slurm.sh
#
set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
RM_INDIR="$STEP02_OUT/repeatmasker"                   # holds the RM_<genome_id>/ dirs (run_repeatmasker.slurm.sh's output)
DROP_LIST="$STEP01_OUT/drop_contigs.txt"               # the step-01 manifest
OUTDIR="$STEP02_OUT/repeatmasker_contig_filtered"     # cleaned RM_<genome_id>/ go here (new dir)
#
# Alternative: instead of the manifest, compare against the final assemblies
# (more accurate .tbl 'total length'). Leave DROP_LIST="" and set:
GENOMES=""                                              # e.g. "genomes"
# ────────────────────────────────────────────────────────────────────────────

find_script PY species_specific/Afum/step02_te_annotation_clean_up/fix_repeatmasker_after_contig_drop.py

### ADAPT ### — put a python3 (>= 3.6, stdlib only) on the PATH. Any env works;
#              the step-01 / synteny env is fine, or the system python3:
# set +u
# module load python/3.11 || true
# set -u
### END ADAPT ###

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
compgen -G "$RM_INDIR"/RM_*/ >/dev/null || compgen -G "$RM_INDIR"/R_*/ >/dev/null \
    || compgen -G "$RM_INDIR"/*.out >/dev/null \
    || { echo "ERROR: no RM_<genome_id>/ folders or *.out files in $RM_INDIR" >&2; exit 1; }

SRC_ARGS=()
if [[ -n "${DROP_LIST:-}" ]]; then
    [[ -f "$DROP_LIST" ]] || { echo "ERROR: drop list not found: $DROP_LIST" >&2; exit 1; }
    SRC_ARGS=(--drop-list "$DROP_LIST")
elif [[ -n "${GENOMES:-}" ]]; then
    [[ -d "$GENOMES" ]] || { echo "ERROR: genomes folder not found: $GENOMES" >&2; exit 1; }
    SRC_ARGS=(--genomes-dir "$GENOMES")
else
    echo "ERROR: set DROP_LIST or GENOMES at the top of the script" >&2; exit 1
fi

python3 "$PY" --rm-indir "$RM_INDIR" "${SRC_ARGS[@]}" --outdir "$OUTDIR"

echo "Done — cleaned RM folders in $OUTDIR/"
