#!/bin/bash
#SBATCH --job-name=get_pangenes_basket
#SBATCH --partition=low
#SBATCH --output=pangenes_basket_%j.log
#SBATCH --error=pangenes_basket_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#
# The same job as run_get_pangenes_job.sh, on the basket annotation (step 05), one mode at a time.
# Needs link_for_get_pangenes_basket.slurm.sh to have already run for that mode. If the output folder
# (<INPUT_DIR>_pangenes) already exists, get_pangenes re-uses the files cached there: rename it before a new run.
# Submit from your project folder, once per mode:
#     sbatch --export=ALL,MODE=keep_similar 06_GET_PANGENES/scripts/run_get_pangenes_basket_job.sh
#     sbatch --export=ALL,MODE=prefer_panel_model 06_GET_PANGENES/scripts/run_get_pangenes_basket_job.sh

set -euo pipefail

# Edit these for your run.
MODE="${MODE:-keep_similar}"                                                    # the basket annotation mode (set with --export)
GET_PANGENES_ENV="/home/aludwig/hc-storage/micromamba/envs/get_pangenes/bin"   # bin/ of the get_pangenes micromamba environment
INPUT_DIR="${INPUT_DIR:-./for_get_pangenes_basket_${MODE}}"   # link_for_get_pangenes_basket.slurm.sh's OUT_DIR
N_CPUS=16                           # match --cpus-per-task above
TAXA_THRESHOLD=0                    # get_pangenes' -t: minimum number of genomes a cluster must appear in to
                                     # be reported at all - 0 reports everything, including genome-private ones
# Remove any previously loaded cluster modules.
module purge

# Avoid mixing Perl libraries from other environments.
unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

# Add the micromamba environment's own programs (get_pangenes.pl, its bundled Python, bedtools) to PATH.
export PATH="$GET_PANGENES_ENV:$PATH"

GETP="$GET_PANGENES_ENV/get_pangenes.pl"
BEDTOOLS="$GET_PANGENES_ENV/bedtools"

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: $INPUT_DIR not found - run link_for_get_pangenes_basket.slurm.sh first" >&2; exit 1; }

if [[ -n "${TAXON_REGEX:-}" ]]; then
  "$GETP" -d "$INPUT_DIR" -t "$TAXA_THRESHOLD" -s "$TAXON_REGEX" -n "$N_CPUS" -B "$BEDTOOLS"
else
  "$GETP" -d "$INPUT_DIR" -t "$TAXA_THRESHOLD" -n "$N_CPUS" -B "$BEDTOOLS"
fi

echo "get_pangenes completed (mode $MODE). Outputs (pangene_matrix*.tab, etc.) are in ${INPUT_DIR}_pangenes"
