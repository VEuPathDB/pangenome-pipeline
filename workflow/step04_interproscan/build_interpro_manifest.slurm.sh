#!/bin/bash
#SBATCH --job-name=build_interpro_manifest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Builds interpro_manifest.tsv: one row per (genome, chunk-of-CHUNK_SIZE-proteins) -
# interproscan_array_job.slurm.sh reads this to know which genome/sequence-range each
# SLURM_ARRAY_TASK_ID should process. Run this FIRST (after extract_proteins.slurm.sh),
# then check the printed "Total chunks: N" line before submitting the array job.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/workflow/step04_interproscan/build_interpro_manifest.slurm.sh
#
set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── ADAPT: the paths below and the #SBATCH lines ───
PROTEIN_DIR="$STEP04_OUT/proteins"    # <genome_label>.faa, from extract_proteins.slurm.sh
ORDER="genome_order.txt"           # step-01 genome list
CHUNK_SIZE=1000
MANIFEST_OUT="$STEP04_OUT/interpro_manifest.tsv"
# ─────────────────────────────────────────────────────

[[ -d "$PROTEIN_DIR" ]] || { echo "ERROR: protein folder not found: $PROTEIN_DIR (run extract_proteins.slurm.sh first)" >&2; exit 1; }
[[ -f "$ORDER" ]] || { echo "ERROR: genome list not found: $ORDER" >&2; exit 1; }
mkdir -p "$(dirname "$MANIFEST_OUT")"

mapfile -t GENOMES < <(grep -vE '^\s*(#|$)' "$ORDER")
[[ ${#GENOMES[@]} -ge 1 ]] || { echo "ERROR: $ORDER has no genomes" >&2; exit 1; }

echo -e "task_id\tgenome\tchunk_num\tstart_seq\tend_seq" > "$MANIFEST_OUT"

task_id=0
for genome in "${GENOMES[@]}"; do
  fasta="$PROTEIN_DIR/${genome}.faa"
  [[ -f "$fasta" ]] || { echo "ERROR: missing $fasta" >&2; exit 1; }

  n_seqs=$(grep -c '^>' "$fasta")
  n_chunks=$(( (n_seqs + CHUNK_SIZE - 1) / CHUNK_SIZE ))
  echo "[$genome] $n_seqs proteins -> $n_chunks chunks"

  for ((c=1; c<=n_chunks; c++)); do
    task_id=$((task_id + 1))
    start=$(( (c - 1) * CHUNK_SIZE + 1 ))
    end=$(( c * CHUNK_SIZE ))
    echo -e "${task_id}\t${genome}\t${c}\t${start}\t${end}" >> "$MANIFEST_OUT"
  done
done

echo ""
echo "Manifest written to $MANIFEST_OUT"
echo "Total chunks: $task_id"
echo "Next: submit interproscan_array_job.slurm.sh with --array=1-${task_id}%<max_concurrent>"
echo "on the sbatch command line (see that script's own header for the exact command) -"
echo "no need to hand-edit its #SBATCH --array line."
