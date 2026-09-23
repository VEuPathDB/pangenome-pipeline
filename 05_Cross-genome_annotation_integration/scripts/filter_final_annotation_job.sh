#!/bin/bash
#SBATCH --job-name=filter_final_annotation
#SBATCH --partition=low
#SBATCH --output=filtered_annotation/logs/%x_%j.out
#SBATCH --error=filtered_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

# Submit from your project folder. The log folder must exist BEFORE submitting (SLURM does not
# create it - without it the job leaves no log):
#     mkdir -p filtered_annotation/logs && sbatch 05_Cross-genome_annotation_integration/scripts/filter_final_annotation_job.sh

# Edit these for your run.
STRAINS="genome_order.txt"                 # comma-separated names or a file with one per line, first = the reference (same as ORDER)
FINAL_DIR="expanded_annotation"            # build_final_annotation_basic.sh's OUTPUT_DIR - the raw, unfiltered <genome>_final_annotation.gff files
OUT_DIR="filtered_annotation"              # where to write the filtered GFFs + reports
MIN_RESCUED_COVERAGE="0.5"                 # every rescued_<genome> gene below this Liftoff coverage is dropped outright
MIN_OVERLAP_BP="20"                        # CDS-to-CDS overlaps under this many bp are ignored - boundary-rounding noise
ID_PREFIX="Afu"                            # renamed gene IDs look like <ID_PREFIX><N>g<00001, ...>

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (standard library only -
### change this only if your cluster needs a specific environment to get a plain python3 on PATH):
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate bioinfo
### END ADAPT ###
set -u

# Locate the .py. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there - fall back to the standard layout
# relative to the directory you submitted from.
PY=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/filter_final_annotation.py" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Cross-genome_annotation_integration/scripts/filter_final_annotation.py" \
    "05_Cross-genome_annotation_integration/scripts/filter_final_annotation.py"
do
    [[ -f "$cand" ]] && { PY="$cand"; break; }
done
[[ -n "$PY" ]] || { echo "ERROR: cannot find filter_final_annotation.py (run from the project folder)" >&2; exit 1; }
mkdir -p "$OUT_DIR"

python3 "$PY" \
  --strains "$STRAINS" --final-dir "$FINAL_DIR" \
  --out-dir "$OUT_DIR" \
  --min-rescued-coverage "$MIN_RESCUED_COVERAGE" \
  --min-overlap-bp "$MIN_OVERLAP_BP" \
  --id-prefix "$ID_PREFIX"

echo "All genomes filtered. Outputs are in $OUT_DIR"
