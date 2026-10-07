#!/bin/bash
#SBATCH --job-name=prepare_for_liftoff
#SBATCH --partition=low
#SBATCH --output=prepare_for_liftoff_%j.out
#SBATCH --error=prepare_for_liftoff_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#
# Prepares each genome's GFF for Liftoff: unifies gene feature types (so Liftoff transfers the
# protein-coding genes and not the non-coding RNA genes) and removes UTRs. See prepare_for_liftoff.py.
#
# Submit from your project folder:
#     sbatch 05_Cross-genome_annotation_integration/scripts/prepare_for_liftoff_job.sh

set -euo pipefail

# ─── Edit these for your run ───
GENOMES="genome_order.txt"                         # comma-separated genome names, or a file with one per line
GFF_DIR="gff_curated"                              # GFFs to prepare (never edited in place) - use gff_with_interpro if no genes were kill-listed
OUT_DIR="gff_for_liftoff"                          # where the prepared GFFs (+ liftoff_feature_types.txt) go
# ────────────────────────────────

set +u
module purge
### ADAPT ### — put python3 on the PATH:
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###
set -u

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script PY 05_Cross-genome_annotation_integration/scripts/prepare_for_liftoff.py


python3 "$PY" \
  --genomes "$GENOMES" --gff-dir "$GFF_DIR" \
  --out-dir "$OUT_DIR"

echo "Done. Prepared GFFs are in $OUT_DIR"
