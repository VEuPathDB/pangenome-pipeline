#!/usr/bin/env bash
#SBATCH --job-name=fastqc
#SBATCH --partition=low
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=fastqc_%j.out
#SBATCH --error=fastqc_%j.err

set -euo pipefail

module load java/17.0.10+7
module load fastqc/0.11.9
module load multiqc/1.14

mkdir -p fastqc_combined

fastqc \
    -t 8 \
    -o fastqc_combined \
    combined_left_trinity.fastq.gz \
    combined_right_trinity.fastq.gz

multiqc fastqc_combined -o multiqc_combined

echo "FastQC and MultiQC completed successfully."
