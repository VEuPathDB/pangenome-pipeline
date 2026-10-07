#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err

set -euo pipefail

# Step 03 jobs are submitted from inside 03_denovo_Genome_Annotation/ (see its README), so look
# for the pipeline's lib/ upwards from there; PIPELINE_DIR overrides.
_pd="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: lib/slurm_common.sh not found above ${SLURM_SUBMIT_DIR:-$PWD} - set PIPELINE_DIR" >&2; exit 1; }
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
