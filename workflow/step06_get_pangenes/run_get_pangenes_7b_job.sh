#!/bin/bash
#SBATCH --job-name=get_pangenes_7b
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=200:00:00
#
# 7b variant (no pre-filter/tag) - see run_get_pangenes_7a_job.sh for the 7a (pre-filter) one.
# Needs link_for_get_pangenes_7b.slurm.sh to have already run.

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
INPUT_DIR="$STEP06_OUT/for_get_pangenes_7b"   # link_for_get_pangenes_7b.slurm.sh's OUT_DIR; results go to <INPUT_DIR>_pangenes/
N_CPUS=32                           # match --cpus-per-task above
TAXA_THRESHOLD=0                    # get_pangenes' -t: minimum number of genomes a cluster must appear in to
                                     # be reported at all - 0 reports everything, including genome-private ones

find_script COMMON workflow/step06_get_pangenes/get_pangenes_common.sh
# shellcheck source=get_pangenes_common.sh
source "$COMMON"

VARIANT="7b"
run_get_pangenes
