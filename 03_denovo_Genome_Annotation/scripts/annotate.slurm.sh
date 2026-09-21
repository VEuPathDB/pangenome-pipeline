#!/usr/bin/env bash
#SBATCH --job-name=fun_annotate
#SBATCH --partition=hi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=annotate_%j.out
#SBATCH --error=annotate_%j.err

set -euo pipefail

SHARED="${FUNANNOTATE_SHARED:-/mnt/hc-storage/groups/veupathdb_annotation/eb_tools/funannotate}"
# Use output from update (if updated) or predict:
GFF_INPUT="funannotate_update/update_results/GENOME.gff3"
FASTA_INPUT="funannotate_update/update_results/GENOME.scaffolds.fa"
SPECIES="Aspergillus fumigatus"
STRAIN="Af293"
OUTDIR="funannotate_annotate"

echo "Starting funannotate annotate: $(date)"

singularity exec \
    -B "${PWD}:${PWD}" \
    -B "${SHARED}/databases/funannotate_db:${SHARED}/databases/funannotate_db" \
    -B "${SHARED}/overlay/usr/share/perl5:/usr/share/perl5" \
    -B "${SHARED}/overlay/usr/local/bin:/usr/local/bin" \
    -B "${SHARED}/overlay/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu" \
    -B "${SHARED}/overlay/venv/lib/python3.8/site-packages/funannotate/config:/venv/lib/python3.8/site-packages/funannotate/config" \
    -B "${SHARED}/pasa.config.txt:/venv/lib/python3.8/site-packages/funannotate/config/pasa.config.txt" \
    --pwd "${PWD}" \
    --env FUNANNOTATE_DB="${SHARED}/databases/funannotate_db" \
    --env PASACONF="/venv/lib/python3.8/site-packages/funannotate/config/pasa.config.txt" \
    --env PATH="/usr/local/bin:/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    "${SHARED}/containers/funannotate.sif" \
    funannotate annotate \
        --gff "${GFF_INPUT}" \
        --fasta "${FASTA_INPUT}" \
        --species "${SPECIES}" \
        --strain "${STRAIN}" \
        --cpus 48 \
        --out "${OUTDIR}"

echo "Finished funannotate annotate at $(date)"
