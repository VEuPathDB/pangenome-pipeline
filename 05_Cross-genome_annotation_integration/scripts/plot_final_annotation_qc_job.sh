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
source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
### ADAPT ### — put micromamba on the PATH, then activate the env (needs Python 3 + matplotlib):
activate_env "$ENV_SYNTENY"  # env names / micromamba: pipeline_config.sh
### END ADAPT ###
set -u

find_script PY 05_Cross-genome_annotation_integration/scripts/plot_final_annotation_qc.py

python3 "$PY" \
  --strains "$STRAINS" --final-dir "$FINAL_DIR" \
  --genome-dir "$GENOME_DIR" --genome-suffix "$GENOME_SUFFIX" \
  --out-dir "$OUT_DIR"

echo "QC plots completed. Outputs are in $OUT_DIR"
