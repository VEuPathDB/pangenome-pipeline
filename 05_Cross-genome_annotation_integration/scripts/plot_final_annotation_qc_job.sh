#!/bin/bash
#SBATCH --job-name=plot_final_annotation_qc
#SBATCH --partition=low
#SBATCH --output=final_annotation_qc_plots/logs/%x_%j.out
#SBATCH --error=final_annotation_qc_plots/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

set -euo pipefail

# Submit from your project folder. The log folder must exist BEFORE submitting (SLURM does not
# create it - without it the job leaves no log):
#     mkdir -p final_annotation_qc_plots/logs && sbatch 05_Cross-genome_annotation_integration/scripts/plot_final_annotation_qc_job.sh

# Edit these for your run.
STRAINS="genome_order.txt"                   # comma-separated names or a file with one per line, first = the reference (same as ORDER)
FINAL_DIR="expanded_annotation"              # build_final_annotation_basic.sh's OUTPUT_DIR
GENOME_DIR="gff_for_liftoff"                 # each genome's own annotation, to compare the gene counts with
GENOME_SUFFIX=".gff3"                        # filename suffix for each genome's own GFF
OUT_DIR="final_annotation_qc_plots"          # where to write the plots and summary.tsv

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (needs Python 3 + matplotlib):
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

# Locate the .py. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there - fall back to the standard layout
# relative to the directory you submitted from.
PY=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/plot_final_annotation_qc.py" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Cross-genome_annotation_integration/scripts/plot_final_annotation_qc.py" \
    "05_Cross-genome_annotation_integration/scripts/plot_final_annotation_qc.py"
do
    [[ -f "$cand" ]] && { PY="$cand"; break; }
done
[[ -n "$PY" ]] || { echo "ERROR: cannot find plot_final_annotation_qc.py (run from the project folder)" >&2; exit 1; }

python3 "$PY" \
  --strains "$STRAINS" --final-dir "$FINAL_DIR" \
  --genome-dir "$GENOME_DIR" --genome-suffix "$GENOME_SUFFIX" \
  --out-dir "$OUT_DIR"

echo "QC plots completed. Outputs are in $OUT_DIR"
