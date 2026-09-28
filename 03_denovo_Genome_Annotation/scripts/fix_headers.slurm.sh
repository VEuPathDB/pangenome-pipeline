#!/usr/bin/env bash
#SBATCH --job-name=fix_headers
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=fix_headers_%j.out
#SBATCH --error=fix_headers_%j.err

set -euo pipefail

module load python/3.11.9
# Step 03 jobs are submitted from inside 03_denovo_Genome_Annotation/ (see its README), so look
# for the pipeline's lib/ upwards from there; PIPELINE_DIR overrides.
_pd="${PIPELINE_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: lib/slurm_common.sh not found above ${SLURM_SUBMIT_DIR:-$PWD} - set PIPELINE_DIR" >&2; exit 1; }

python3 -m pangenome_pipeline.step03_denovo_genome_annotation.fix_headers \
    combined_left_trinity.fastq.gz \
    combined_right_trinity.fastq.gz \
    combined_left_trinity.fixed_headers.fastq.gz \
    combined_right_trinity.fixed_headers.fastq.gz

echo "Headers fixed successfully."
