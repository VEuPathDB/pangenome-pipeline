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

python3 scripts/fix_headers.py \
    combined_left_trinity.fastq.gz \
    combined_right_trinity.fastq.gz \
    combined_left_trinity.fixed_headers.fastq.gz \
    combined_right_trinity.fixed_headers.fastq.gz

echo "Headers fixed successfully."
