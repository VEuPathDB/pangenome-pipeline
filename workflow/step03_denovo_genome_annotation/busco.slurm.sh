#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
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
require_config FUNANNOTATE_SHARED
SHARED="$FUNANNOTATE_SHARED"   # shared funannotate install, set in pipeline_config.sh
PROTEINS_FASTA="funannotate_annotate/annotate_results/GENOME.proteins.fa"
LINEAGE="${SHARED}/databases/funannotate_db/dikarya" # or ascomycota / fungi_odb10
OUTNAME="busco_assessment"

module purge
module load busco/5.4.6

echo "Starting BUSCO evaluation: $(date)"

busco \
    -i "${PROTEINS_FASTA}" \
    -m proteins \
    -l "${LINEAGE}" \
    -o "${OUTNAME}" \
    -c 32

echo "Finished BUSCO evaluation at $(date)"
