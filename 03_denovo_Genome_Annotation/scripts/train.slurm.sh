#!/usr/bin/env bash
#SBATCH --job-name=fun_train
#SBATCH --partition=hi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=train_%j.out
#SBATCH --error=train_%j.err

set -euo pipefail

# Step 03 jobs are submitted from inside 03_denovo_Genome_Annotation/ (see its README), so look
# for the pipeline's lib/ upwards from there; PIPELINE_DIR overrides.
_pd="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: lib/slurm_common.sh not found above ${SLURM_SUBMIT_DIR:-$PWD} - set PIPELINE_DIR" >&2; exit 1; }
require_config FUNANNOTATE_SHARED
SHARED="$FUNANNOTATE_SHARED"   # shared funannotate install, set in pipeline_config.sh
GENOME_MASKED="funannotate/repeatmask/GENOME.masked.fna"
LEFT_READS="combined_left_trinity.fixed_headers.fastq.gz"
RIGHT_READS="combined_right_trinity.fixed_headers.fastq.gz"
SPECIES="${SPECIES:?set SPECIES here (e.g. \"Aspergillus fumigatus\") or with sbatch --export=ALL,SPECIES=...}"
STRAIN="${STRAIN:?set STRAIN here (e.g. Af293) or with sbatch --export=ALL,STRAIN=...}"
OUTDIR="funannotate_train"

export _JAVA_OPTIONS="-Xmx32g"

echo "Starting funannotate train: $(date)"

singularity exec \
    -B "${PWD}:${PWD}" \
    -B "${SHARED}/databases/funannotate_db:${SHARED}/databases/funannotate_db" \
    -B "${SHARED}/overlay/usr/share/perl5:/usr/share/perl5" \
    -B "${SHARED}/overlay/usr/local/bin:/usr/local/bin" \
    -B "${SHARED}/overlay/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu" \
    -B "${SHARED}/pasa.config.txt:/venv/lib/python3.8/site-packages/funannotate/config/pasa.config.txt" \
    --pwd "${PWD}" \
    --env FUNANNOTATE_DB="${SHARED}/databases/funannotate_db" \
    --env PASACONF="/venv/lib/python3.8/site-packages/funannotate/config/pasa.config.txt" \
    --env GENEMARK_PATH="/usr/local/bin/genemark" \
    --env PATH="/usr/local/bin:/venv/bin:/usr/bin:/bin" \
    --env _JAVA_OPTIONS="-Xmx32g" \
    "${SHARED}/containers/funannotate.sif" \
    funannotate train \
        -i "${GENOME_MASKED}" \
        -o "${OUTDIR}" \
        --left "${LEFT_READS}" \
        --right "${RIGHT_READS}" \
        --species "${SPECIES}" \
        --strain "${STRAIN}" \
        --cpus 48 \
        --max_intronlen 3000

echo "Finished funannotate train at $(date)"
