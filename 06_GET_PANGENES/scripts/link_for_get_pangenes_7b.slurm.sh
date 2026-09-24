#!/bin/bash
#SBATCH --job-name=link_for_get_pangenes_7b
#SBATCH --partition=low
#SBATCH --output=link_for_get_pangenes_7b_%j.out
#SBATCH --error=link_for_get_pangenes_7b_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#
# 7b variant (no pre-filter/tag) - see link_for_get_pangenes_7a.slurm.sh for the 7a (pre-filter)
# one. Both can be submitted and run independently; each writes to its own OUT_DIR so a
# get_pangenes run started from one never touches the other's input.
#
# Replaces the old "Step 1: Prepare" - step 05's 7b already writes get_pangenes-ready GFFs
# directly (feature type "gene", IDs unique across the whole panel), so nothing needs
# relabeling any more. All that's actually left before get_pangenes can run: every genome's
# GFF and its matching genome FASTA need to sit in the same directory with the exact same
# basename (get_pangenes.pl pairs them up by filename). This job just symlinks both into
# place - no file is copied or rewritten, so it's instant even for a lot of genomes.
#
# Submit from your project folder:
#     sbatch 06_GET_PANGENES/scripts/link_for_get_pangenes_7b.slurm.sh

set -euo pipefail

# Edit these for your run.
ANNOTATION_DIR="tagged_annotation"                    # step 05's 7b OUT_DIR
ANNOTATION_SUFFIX="_final_annotation_tagged.gff"       # matches 7b's naming
GENOME_DIR="genomes"                                   # step-01 assemblies folder
GENOME_SUFFIX=".fna"                                   # filename suffix for each genome's own FASTA
ORDER="genome_order.txt"                               # same genome list used throughout this pipeline
OUT_DIR="for_get_pangenes_7b"                          # point run_get_pangenes_7b_job.sh's INPUT_DIR here

if [[ -f "$ORDER" ]]; then
  mapfile -t GENOMES < <(grep -vE '^\s*(#|$)' "$ORDER")
else
  IFS=',' read -ra GENOMES <<< "$ORDER"
fi
[[ ${#GENOMES[@]} -ge 1 ]] || { echo "ERROR: ORDER has no genomes - check it's been edited" >&2; exit 1; }

mkdir -p "$OUT_DIR"

n=0
for genome in "${GENOMES[@]}"; do
  gff="$ANNOTATION_DIR/${genome}${ANNOTATION_SUFFIX}"
  fna="$GENOME_DIR/${genome}${GENOME_SUFFIX}"
  [[ -f "$gff" ]] || { echo "ERROR: missing $gff" >&2; exit 1; }
  [[ -f "$fna" ]] || { echo "ERROR: missing $fna" >&2; exit 1; }
  ln -sf "$(cd "$(dirname "$gff")" && pwd)/$(basename "$gff")" "$OUT_DIR/${genome}_for_get_pangenes.gff"
  ln -sf "$(cd "$(dirname "$fna")" && pwd)/$(basename "$fna")" "$OUT_DIR/${genome}_for_get_pangenes.fna"
  n=$((n + 1))
  echo "linked $genome"
done

echo "$n genome(s) linked into $OUT_DIR. Outputs are in $OUT_DIR"
