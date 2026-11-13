#!/bin/bash
#SBATCH --job-name=link_for_get_pangenes_basket
#SBATCH --partition=low
#SBATCH --output=link_for_get_pangenes_basket_%j.out
#SBATCH --error=link_for_get_pangenes_basket_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#
# The same job as link_for_get_pangenes.slurm.sh, for the GFFs of the basket annotation (step 05,
# 05_Basket_annotation_integration), one mode at a time. Those GFFs are get_pangenes-ready (feature type "gene",
# IDs unique across the panel), but they hold every contig of each genome's own annotation: a gene on a contig
# that is not in the genome FASTA used here (e.g. a removed redundant contig) is left out, and counted in the log.
# The genome FASTA is symlinked; the GFF is copied without those genes.
#
# Submit from your project folder (MODE=keep_similar for the other step 05 mode, if it was built):
#     sbatch 06_GET_PANGENES/scripts/link_for_get_pangenes_basket.slurm.sh

set -euo pipefail

# Edit these for your run.
MODE="${MODE:-prefer_panel_model}"                                  # the basket annotation mode (set with --export)
ANNOTATION_DIR="${ANNOTATION_DIR:-basket_annotation/gff_${MODE}_final}"   # the final GFFs of step 05 (after miniprot and the overlap rule)
ANNOTATION_SUFFIX="_basket_annotation.gff"
GENOME_DIR="genomes"                                          # step-01 assemblies folder
GENOME_SUFFIX=".fna"                                          # filename suffix for each genome's own FASTA
ORDER="genome_order.txt"                                      # same genome list used throughout this pipeline
OUT_DIR="${OUT_DIR:-for_get_pangenes_basket_${MODE}}"         # run_get_pangenes_basket_job.sh's INPUT_DIR

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
  out="$OUT_DIR/${genome}_for_get_pangenes.gff"
  # keep only the lines on contigs of the FASTA (a gene and its mRNA/exons/CDS share the contig)
  awk -F'\t' 'NR == FNR { if (/^>/) { sub(/^>/, ""); split($0, a, /[ \t]/); keep[a[1]] = 1 } next }
              /^#/ || ($1 in keep) { print; next }
              $3 == "gene" { dropped++ }
              END { print dropped + 0 > "/dev/stderr" }' "$fna" "$gff" > "$out" 2> "$OUT_DIR/.dropped"
  ln -sf "$(cd "$(dirname "$fna")" && pwd)/$(basename "$fna")" "$OUT_DIR/${genome}_for_get_pangenes.fna"
  n=$((n + 1))
  echo "$genome: $(grep -c $'\tgene\t' "$out") genes; left out (contig not in the FASTA): $(cat "$OUT_DIR/.dropped")"
done
rm -f "$OUT_DIR/.dropped"

echo "$n genome(s) prepared in $OUT_DIR"
