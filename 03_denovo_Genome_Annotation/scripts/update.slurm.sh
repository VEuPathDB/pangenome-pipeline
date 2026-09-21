#!/usr/bin/env bash
#SBATCH --job-name=fun_update
#SBATCH --partition=hi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=update_%j.out
#SBATCH --error=update_%j.err

set -euo pipefail

SHARED="${FUNANNOTATE_SHARED:-/mnt/hc-storage/groups/veupathdb_annotation/eb_tools/funannotate}"
GENOME_FASTA="funannotate_predict/predict_results/GENOME.scaffolds.fa"
GBK_INPUT="funannotate_annotate/annotate_results/GENOME.gbk"
# Or from predict if skipping initial annotate:
# GBK_INPUT="funannotate_predict/predict_results/GENOME.gbk"
TRINITY_FASTA="funannotate_train/training/trinity.fasta.clean"
LEFT_NORM="funannotate_train/training/left.fq.gz"
RIGHT_NORM="funannotate_train/training/right.fq.gz"
OUTDIR="funannotate_update"

export _JAVA_OPTIONS="-Xmx32g"

echo "Starting funannotate update: $(date)"

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
    --env _JAVA_OPTIONS="-Xmx32g" \
    "${SHARED}/containers/funannotate.sif" \
    funannotate update \
        --input "${GBK_INPUT}" \
        --fasta "${GENOME_FASTA}" \
        --trinity "${TRINITY_FASTA}" \
        --left_norm "${LEFT_NORM}" \
        --right_norm "${RIGHT_NORM}" \
        --no_trimmomatic \
        --no_normalize_reads \
        --out "${OUTDIR}" \
        --cpus 48 \
        --max_intronlen 3000

echo "Finished funannotate update at $(date)"
