#!/bin/bash
#SBATCH --job-name=check_tandem_pairs
#SBATCH --partition=low
#SBATCH --output=basket_annotation/logs/%x_%j.out
#SBATCH --error=basket_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

# Check of the tandem-array detection, after basket_annotation_job.sh: every neighbouring pair compared by
# baskets_v2.py, with its blastp identity and coverage, and an independent blastp of the closest pairs.
# Writes basket_annotation/tandem_pairs_check.tsv; the summary is in the .out log.
# Submit from your project folder:
#     sbatch 05_Cross-genome_annotation_integration/scripts/check_tandem_pairs_job.sh

# ===== Edit these (the same values as in basket_annotation_job.sh) =====
GENOMES="${GENOMES-genome_order.txt}"
OUT_DIR="${OUT_DIR-basket_annotation}"
GENOME_DIR="${GENOME_DIR-genomes}"
GENOME_SUFFIX="${GENOME_SUFFIX-.fna}"
THREADS="${THREADS-4}"                                          # match --cpus-per-task above
# =======================================================================

set +u
module purge
### ADAPT ### — put python3, blastp and makeblastdb on the PATH:
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
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Cross-genome_annotation_integration/scripts" \
    "05_Cross-genome_annotation_integration/scripts"
do
    [[ -f "$cand/check_tandem_pairs.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find check_tandem_pairs.py (run from the project folder)" >&2; exit 1; }
[[ -d "$OUT_DIR/input/native" ]] || { echo "ERROR: $OUT_DIR/input/native not found: run basket_annotation_job.sh first" >&2; exit 1; }

python3 "$SD/check_tandem_pairs.py" --baskets "$OUT_DIR/baskets_v2.tsv" --native-dir "$OUT_DIR/input/native" \
  --genomes "$GENOMES" --genome-dir "$GENOME_DIR" --genome-suffix "$GENOME_SUFFIX" --threads "$THREADS" \
  --out "$OUT_DIR/tandem_pairs_check"
