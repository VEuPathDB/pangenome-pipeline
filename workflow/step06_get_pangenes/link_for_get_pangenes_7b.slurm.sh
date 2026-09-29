#!/bin/bash
#SBATCH --job-name=link_for_get_pangenes_7b
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#
# 7b variant (no pre-filter/tag) - see link_for_get_pangenes_7a.slurm.sh for the 7a (pre-filter)
# one. Both can be submitted and run independently; each writes to its own OUT_DIR so a
# get_pangenes run started from one never touches the other's input.
#
# Replaces the old "Step 1: Prepare" - step 05's 7b already writes get_pangenes-ready GFFs
# directly (feature type "gene", IDs unique across the whole panel), so nothing needs
# relabeling any more. All that's actually left before get_pangenes can run: every genome's
# GFF and its matching genome FASTA need to sit in the same directory with the exact same
# basename (get_pangenes.pl pairs them up by filename). This job just symlinks both into
# place - no file is copied or rewritten, so it's instant even for a lot of genomes.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/workflow/step06_get_pangenes/link_for_get_pangenes_7b.slurm.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
ANNOTATION_DIR="$STEP05_OUT/tagged_annotation"                    # step 05's 7b OUT_DIR
ANNOTATION_SUFFIX="_final_annotation_tagged.gff"       # matches 7b's naming
GENOME_DIR="genomes"                                   # step-01 assemblies folder
GENOME_SUFFIX=".fna"                                   # filename suffix for each genome's own FASTA
ORDER="genome_order.txt"                               # same genome list used throughout this pipeline
OUT_DIR="$STEP06_OUT/for_get_pangenes_7b"                          # point run_get_pangenes_7b_job.sh's INPUT_DIR here

find_script COMMON workflow/step06_get_pangenes/get_pangenes_common.sh
# shellcheck source=get_pangenes_common.sh
source "$COMMON"

link_for_get_pangenes
