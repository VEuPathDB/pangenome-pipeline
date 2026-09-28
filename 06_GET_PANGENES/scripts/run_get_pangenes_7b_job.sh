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

# Find the shared body. Under sbatch the job script is copied to a spool dir, so
# ${BASH_SOURCE[0]} isn't next to it - fall back to the layout under the submit dir.
COMMON=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/get_pangenes_common.sh" \
    "${SLURM_SUBMIT_DIR:-$PWD}/06_GET_PANGENES/scripts/get_pangenes_common.sh"
do
  [[ -f "$cand" ]] && { COMMON="$cand"; break; }
done
[[ -n "$COMMON" ]] || { echo "ERROR: cannot find get_pangenes_common.sh (submit from the project folder)" >&2; exit 1; }
# shellcheck source=get_pangenes_common.sh
source "$COMMON"

VARIANT="7b"
run_get_pangenes
