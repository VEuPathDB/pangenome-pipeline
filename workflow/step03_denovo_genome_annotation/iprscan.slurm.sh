#!/usr/bin/env bash
#SBATCH --job-name=iprscan
#SBATCH --partition=hi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }
# Step 03 works in one folder, $STEP03_OUT: the inputs you put there (genomes/genome.fna, the
# RNA-seq reads, curated_TE_library.fa) and every output. The paths below are relative to it.
mkdir -p "$STEP03_OUT" && cd "$STEP03_OUT"

PROTEINS_FASTA="funannotate_annotate/annotate_results/GENOME.proteins.fa"
OUTPUT_XML="iprscan.xml"

module purge
module load adoptopenjdk/11.0.12+7
module load interproscan/5.66-98.0

echo "Starting InterProScan: $(date)"

interproscan.sh \
    -i "${PROTEINS_FASTA}" \
    -f XML \
    -o "${OUTPUT_XML}" \
    --goterms \
    --iprlookup \
    --pathways \
    -cpu 48

echo "Finished InterProScan at $(date)"
