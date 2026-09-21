#!/usr/bin/env bash
#SBATCH --job-name=iprscan
#SBATCH --partition=hi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=iprscan_%j.out
#SBATCH --error=iprscan_%j.err

set -euo pipefail

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
