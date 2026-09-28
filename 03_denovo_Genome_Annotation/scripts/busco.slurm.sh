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

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # pipeline_config.sh
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
