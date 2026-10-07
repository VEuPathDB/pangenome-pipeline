#!/bin/bash
#SBATCH --job-name=get_pangenes_7a
#SBATCH --partition=low
#SBATCH --output=pangenes_7a_%j.log
#SBATCH --error=pangenes_7a_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=200:00:00
#
# 7a variant (pre-filter) - see run_get_pangenes_7b_job.sh for the 7b (no pre-filter/tag) one.
# Needs link_for_get_pangenes_7a.slurm.sh to have already run.

set -euo pipefail

# Edit these for your run.
INPUT_DIR="./for_get_pangenes_7a"   # link_for_get_pangenes_7a.slurm.sh's OUT_DIR
N_CPUS=32                           # match --cpus-per-task above
TAXA_THRESHOLD=0                    # get_pangenes' -t: minimum number of genomes a cluster must appear in to
                                     # be reported at all - 0 reports everything, including genome-private ones

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script COMMON 06_GET_PANGENES/scripts/get_pangenes_common.sh
# shellcheck source=get_pangenes_common.sh
source "$COMMON"

VARIANT="7a"
run_get_pangenes
