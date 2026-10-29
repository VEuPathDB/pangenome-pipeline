#!/bin/bash
#SBATCH --job-name=basket_member_psauron
#SBATCH --partition=low
#SBATCH --output=basket_annotation/psauron_members/logs/%x_%j.out
#SBATCH --error=basket_annotation/psauron_members/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Evaluation of the winner rule of step 05 (nothing is changed): the PSAURON score of every member of the baskets with
# two or more winners, and of the baskets whose reference winner has a valid ORF in at most half of the genomes where it
# was placed. Run after basket_annotation_job.sh. The log folder must exist BEFORE submitting:
#     mkdir -p basket_annotation/psauron_members/logs && sbatch 05_Basket_annotation_integration/scripts/basket_member_psauron_job.sh

# ===== Edit these =====
BASKET_DIR="${BASKET_DIR-basket_annotation}"                 # step 05 output
NATIVE_DIR="${NATIVE_DIR-basket_annotation/input/native}"    # the native GFFs step 05 read
PSAURON_DIR="${PSAURON_DIR-psauron}"                         # step 04 PSAURON scores (*.csv)
OUT_DIR="${OUT_DIR-basket_annotation/psauron_members}"
# ======================

set +u
module purge 2>/dev/null || true
### ADAPT ### — put python3 on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

# Locate the scripts. Under `sbatch` the wrapper is copied to a spool dir, so ${BASH_SOURCE[0]} is useless there -
# fall back to the standard layouts relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Basket_annotation_integration/scripts" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Cross-genome_annotation_integration/scripts" \
    "05_Basket_annotation_integration/scripts" \
    "05_Cross-genome_annotation_integration/scripts"
do
    [[ -f "$cand/basket_member_psauron.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find basket_member_psauron.py (run from the project folder)" >&2; exit 1; }
for d in "$BASKET_DIR" "$NATIVE_DIR" "$PSAURON_DIR"; do
    [[ -d "$d" ]] || { echo "ERROR: missing folder $d" >&2; exit 1; }
done

PYTHONPATH="$SD" python3 "$SD/basket_member_psauron.py" --basket-dir "$BASKET_DIR" --native-dir "$NATIVE_DIR" \
  --psauron-dir "$PSAURON_DIR" --out-dir "$OUT_DIR"

echo "Done. $OUT_DIR/basket_psauron_summary.tsv and basket_members_psauron.tsv"
