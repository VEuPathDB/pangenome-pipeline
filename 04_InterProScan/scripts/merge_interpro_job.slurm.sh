#!/bin/bash
#SBATCH --job-name=merge_interpro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=merge_interpro_%j.out
#SBATCH --error=merge_interpro_%j.err
#
# Only run this once EVERY interproscan_array_job.slurm.sh task has finished - a genome
# with a missing chunk will still produce a lookup table, just an incomplete one; the
# completeness check below (merge_interpro_results.py) warns about it rather than fail
# loudly, so check each genome's log line after this runs.
#
# For every genome, collapses its interproscan_results/<genome>/chunk_*_interproscan.tsv
# files into one interpro_lookup/<genome>_interpro.tsv (protein_id -> InterPro/Pfam IDs).
#
# Submit from your project folder:
#     sbatch 04_InterProScan/scripts/merge_interpro_job.slurm.sh
#
# ─── ADAPT: the paths below and the #SBATCH lines ───
RESULTS_DIR="interpro/interproscan_results"
PROTEIN_DIR="interpro/proteins"     # for the completeness check
ORDER="genome_order.txt"
OUT_DIR="interpro/interpro_lookup"
# ─────────────────────────────────────────────────────

set -euo pipefail

# Locate the .py. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there — fall back to the standard layout
# relative to the directory you submitted from.
PY=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/merge_interpro_results.py" \
    "${SLURM_SUBMIT_DIR:-$PWD}/04_InterProScan/scripts/merge_interpro_results.py" \
    "04_InterProScan/scripts/merge_interpro_results.py"
do
    [[ -f "$cand" ]] && { PY="$cand"; break; }
done
[[ -n "$PY" ]] || { echo "ERROR: cannot find merge_interpro_results.py (run from the project folder)" >&2; exit 1; }

[[ -d "$RESULTS_DIR" ]] || { echo "ERROR: results folder not found: $RESULTS_DIR (run the array job first)" >&2; exit 1; }
[[ -f "$ORDER" ]] || { echo "ERROR: genome list not found: $ORDER" >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }

mkdir -p "$OUT_DIR"

mapfile -t GENOMES < <(grep -vE '^\s*(#|$)' "$ORDER")

i=0
for genome in "${GENOMES[@]}"; do
  i=$((i + 1))
  echo "=== Genome $i/${#GENOMES[@]}: $genome ==="
  python3 "$PY" \
      --results-dir "$RESULTS_DIR/${genome}" \
      --out "$OUT_DIR/${genome}_interpro.tsv" \
      --fasta "$PROTEIN_DIR/${genome}.faa"
done

echo
echo "Done — all genomes merged, output in $OUT_DIR/"
