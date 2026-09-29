#!/usr/bin/env bash
#SBATCH --job-name=fun_train
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
