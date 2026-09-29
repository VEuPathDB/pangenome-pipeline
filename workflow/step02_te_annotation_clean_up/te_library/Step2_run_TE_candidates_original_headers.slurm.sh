#!/usr/bin/env bash
#SBATCH --job-name=TE_step2_candidates
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --partition=low

set -euo pipefail

step_start="$(date +%s)"

# Submit from the project folder (mkdir -p logs first):
#   sbatch $PIPELINE_DIR/workflow/step02_te_annotation_clean_up/te_library/Step2_run_TE_candidates_original_headers.slurm.sh
# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }
SLURM_DIR="${PIPELINE_DIR}/workflow/step02_te_annotation_clean_up/te_library"
source "${SLURM_DIR}/config.sh"
cd "${PROJECT_DIR}"

mkdir -p "${LOG_DIR}" "${CANDIDATE_OUTDIR}"

echo "============================================================"
echo "Step 2: remove known TE families"
echo "Started: $(date)"
echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "Project directory: ${PROJECT_DIR}"
echo "============================================================"

summary="${CANDIDATE_OUTDIR}/summary.tsv"
if [[ -s "${summary}" ]]; then
    echo "Step 2 output already exists; skipping."
    echo "Existing summary: ${summary}"
else
    load_filtering_environment

    "${SLURM_DIR}/TE_candidate_pipeline_original_headers.sh" \
        --input_dir "${REPEATMODELER_OUTDIR}" \
        --library "${CURATED_TE_LIBRARY}" \
        --output_dir "${CANDIDATE_OUTDIR}" \
        --consensi_name "consensi.fa.classified" \
        --identity "${KNOWN_TE_MIN_IDENTITY}" \
        --coverage "${KNOWN_TE_MIN_QCOV}" \
        --min_length "${KNOWN_TE_MIN_LENGTH}" \
        --threads "${FILTER_THREADS}"
fi

step_end="$(date +%s)"
echo "Step 2 finished: $(date)"
echo "Elapsed seconds: $((step_end - step_start))"
