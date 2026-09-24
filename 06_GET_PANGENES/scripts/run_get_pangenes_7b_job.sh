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

if [[ -n "${TAXON_REGEX:-}" ]]; then
  "$GETP" -d "$INPUT_DIR" -t "$TAXA_THRESHOLD" -s "$TAXON_REGEX" -n "$N_CPUS" -B "$BEDTOOLS"
else
  "$GETP" -d "$INPUT_DIR" -t "$TAXA_THRESHOLD" -n "$N_CPUS" -B "$BEDTOOLS"
fi

echo "get_pangenes (7b) completed. Outputs (pangene_matrix*.tab, etc.) are in $INPUT_DIR"
