#!/bin/bash
#SBATCH --job-name=get_pangenes_7b
#SBATCH --partition=low
#SBATCH --output=pangenes_7b_%j.log
#SBATCH --error=pangenes_7b_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=200:00:00
#
# 7b variant (no pre-filter/tag) - see run_get_pangenes_7a_job.sh for the 7a (pre-filter) one.
# Needs link_for_get_pangenes_7b.slurm.sh to have already run.

set -euo pipefail

# Edit these for your run.
GET_PANGENES_ENV="/home/aludwig/hc-storage/micromamba/envs/get_pangenes/bin"   # bin/ of the get_pangenes micromamba environment
INPUT_DIR="./for_get_pangenes_7b"   # link_for_get_pangenes_7b.slurm.sh's OUT_DIR
N_CPUS=32                           # match --cpus-per-task above
TAXA_THRESHOLD=0                    # get_pangenes' -t: minimum number of genomes a cluster must appear in to
                                     # be reported at all - 0 reports everything, including genome-private ones

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script COMMON 06_GET_PANGENES/scripts/get_pangenes_common.sh
# shellcheck source=get_pangenes_common.sh
source "$COMMON"

VARIANT="7b"
run_get_pangenes
