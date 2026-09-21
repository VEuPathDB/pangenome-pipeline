#!/usr/bin/env bash
#SBATCH --job-name=repeatmask
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=repeatmask_%j.out
#SBATCH --error=repeatmask_%j.err

set -euo pipefail

SHARED="${FUNANNOTATE_SHARED:-/mnt/hc-storage/groups/veupathdb_annotation/eb_tools/funannotate}"
GENOME_IN="funannotate/GENOME.sorted.fna"
TE_LIB="curated_TE_library.fa"
OUTDIR="funannotate/repeatmask"

mkdir -p "${OUTDIR}"

echo "Starting RepeatMasker: $(date)"

# Option A: Using RepeatMasker directly
singularity exec \
    -B "${PWD}:${PWD}" \
    "${SHARED}/containers/dfam-tetools-latest.sif" \
    RepeatMasker \
        -pa 32 \
        -lib "${TE_LIB}" \
        -dir "${OUTDIR}" \
        -xsmall \
        -gff \
        -html \
        "${GENOME_IN}"

# Option B (alternative): Using funannotate mask
# singularity exec \
#     -B "${PWD}:${PWD}" \
#     -B "${SHARED}/databases/funannotate_db:${SHARED}/databases/funannotate_db" \
#     --env FUNANNOTATE_DB="${SHARED}/databases/funannotate_db" \
#     "${SHARED}/containers/funannotate.sif" \
#     funannotate mask \
#         -i "${GENOME_IN}" \
#         -o "${OUTDIR}/GENOME.masked.fna" \
#         --cpus 32

echo "Finished repeat masking at $(date)"
