#!/bin/bash
#SBATCH --job-name=link_for_get_pangenes_7a
#SBATCH --partition=low
#SBATCH --output=link_for_get_pangenes_7a_%j.out
#SBATCH --error=link_for_get_pangenes_7a_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#
# 7a variant (pre-filter) - see link_for_get_pangenes_7b.slurm.sh for the 7b (no pre-filter/tag)
# one. Both can be submitted and run independently; each writes to its own OUT_DIR so a
# get_pangenes run started from one never touches the other's input.
#
# Replaces the old "Step 1: Prepare" - step 05's 7a already writes get_pangenes-ready GFFs
# directly (feature type "gene", IDs unique across the whole panel), so nothing needs
# relabeling any more. All that's actually left before get_pangenes can run: every genome's
# GFF and its matching genome FASTA need to sit in the same directory with the exact same
# basename (get_pangenes.pl pairs them up by filename). This job just symlinks both into
# place - no file is copied or rewritten, so it's instant even for a lot of genomes.
#
# Submit from your project folder:
#     sbatch 06_GET_PANGENES/scripts/link_for_get_pangenes_7a.slurm.sh

set -euo pipefail

# Edit these for your run.
ANNOTATION_DIR="filtered_annotation"                  # step 05's 7a OUT_DIR
ANNOTATION_SUFFIX="_final_annotation_filtered.gff"     # matches 7a's naming
GENOME_DIR="genomes"                                   # step-01 assemblies folder
GENOME_SUFFIX=".fna"                                   # filename suffix for each genome's own FASTA
ORDER="genome_order.txt"                               # same genome list used throughout this pipeline
OUT_DIR="for_get_pangenes_7a"                          # point run_get_pangenes_7a_job.sh's INPUT_DIR here

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script COMMON 06_GET_PANGENES/scripts/get_pangenes_common.sh
# shellcheck source=get_pangenes_common.sh
source "$COMMON"

link_for_get_pangenes
