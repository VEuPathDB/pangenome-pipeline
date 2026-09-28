#!/bin/bash
#SBATCH --job-name=analyze_input_gffs
#SBATCH --partition=low
#SBATCH --output=analyze_input_gffs_%j.out
#SBATCH --error=analyze_input_gffs_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#
# Profiles every genome's own input GFF before Liftoff runs - see analyze_input_gffs.py's own
# header for exactly what it checks and what it writes.
#
# Submit from your project folder:
#     sbatch 05_Cross-genome_annotation_integration/scripts/analyze_input_gffs_job.sh

set -euo pipefail

# ─── Edit these for your run ───
GENOMES="genome_order.txt"   # comma-separated genome names, or a file with one per line
GFF_DIR="gff_with_interpro"   # step 04 output (domain evidence written in) - run this check once, after InterProScan
OUT_DIR="annotation_check"   # where to write the tables and plots
# ────────────────────────────────

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (needs Python 3 + matplotlib):
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script PY 05_Cross-genome_annotation_integration/scripts/analyze_input_gffs.py

python3 "$PY" \
  --genomes "$GENOMES" --gff-dir "$GFF_DIR" --out-dir "$OUT_DIR"

echo "Annotation check completed. Outputs are in $OUT_DIR"
