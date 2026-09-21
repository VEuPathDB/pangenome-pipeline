#!/usr/bin/env bash
#SBATCH --job-name=fun_predict
#SBATCH --partition=hi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=predict_%j.out
#SBATCH --error=predict_%j.err

set -euo pipefail

SHARED="${FUNANNOTATE_SHARED:-/mnt/hc-storage/groups/veupathdb_annotation/eb_tools/funannotate}"
GENOME_MASKED="funannotate/repeatmask/GENOME.masked.fna"
SPECIES="Aspergillus fumigatus"
STRAIN="Af293"
OUTDIR="funannotate_predict"

# Augustus requires a writable config directory
AUG_CONFIG="${PWD}/augustus_config_predict"
mkdir -p "${AUG_CONFIG}"
if [ ! -f "${AUG_CONFIG}/model/states_shadow.cfg" ]; then
    echo "Creating writable Augustus configuration..."
    singularity exec \
        -B "${AUG_CONFIG}":/aug_out \
        "${SHARED}/containers/funannotate.sif" \
        cp -a /usr/share/augustus/config/. /aug_out/
fi

echo "Starting funannotate predict: $(date)"

singularity exec \
    -B "${PWD}:${PWD}" \
    -B "${SHARED}/databases/funannotate_db:${SHARED}/databases/funannotate_db" \
    -B "${SHARED}/overlay/usr/share/perl5:/usr/share/perl5" \
    -B "${SHARED}/overlay/usr/local/bin:/usr/local/bin" \
    -B "${SHARED}/overlay/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu" \
    -B "${SHARED}/overlay/usr/local/bin/genemark:/usr/local/bin/genemark" \
    -B "${SHARED}/pasa.config.txt:/venv/lib/python3.8/site-packages/funannotate/config/pasa.config.txt" \
    -B "${SHARED}/overlay/venv/lib/python3.8/site-packages/funannotate/config:/venv/lib/python3.8/site-packages/funannotate/config" \
    -B "${AUG_CONFIG}:/usr/share/augustus/config" \
    --pwd "${PWD}" \
    --env FUNANNOTATE_DB="${SHARED}/databases/funannotate_db" \
    --env PASACONF="/venv/lib/python3.8/site-packages/funannotate/config/pasa.config.txt" \
    --env GENEMARK_PATH="/usr/local/bin/genemark" \
    --env PATH="/usr/local/bin:/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin" \
    "${SHARED}/containers/funannotate.sif" \
    funannotate predict \
        -i "${GENOME_MASKED}" \
        -o "${OUTDIR}" \
        --species "${SPECIES}" \
        --strain "${STRAIN}" \
        --cpus 48 \
        --organism fungus \
        --max_intronlen 3000

echo "Finished funannotate predict at $(date)"
