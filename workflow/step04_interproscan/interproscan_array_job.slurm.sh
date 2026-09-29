#!/bin/bash
#SBATCH --job-name=interproscan_array
#SBATCH --array=1-1%20    # Smoke-test-only default (1 task) - DO NOT hand-edit this line to
                           # the real chunk count. Always submit with `--array=...` on the
                           # sbatch command line instead (see "Submit from your project
                           # folder" below) - it overrides this line automatically, computed
                           # fresh from the manifest every time, so this can never go stale.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=98:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# Runs InterProScan on one (genome, chunk) pair per array task, as defined by
# interpro_manifest.tsv (built by build_interpro_manifest.slurm.sh - run that first).
# Output: $STEP04_OUT/interproscan_results/<genome>/chunk_<N>_interproscan.tsv
#
# The compute nodes on this cluster can't reach EBI, so -dp (disable precalculated
# match lookup) is required below - every sequence is computed locally instead of
# checked against EBI's online cache. This makes each chunk slower than a
# lookup-enabled run, already accounted for in --time above.
#
# Submit from your project folder - this reads the real chunk count straight out of the
# manifest and passes it to `--array` on the command line, which overrides the smoke-test
# default in the #SBATCH line above, so there is nothing to hand-edit:
#     N=$(($(wc -l < output/step04_interproscan/interpro_manifest.tsv) - 1))   # -1 for the header row
#     # (the manifest is in $STEP04_OUT - by default output/step04_interproscan/ in the project)
#     sbatch --array=1-${N}%20 $PIPELINE_DIR/workflow/step04_interproscan/interproscan_array_job.slurm.sh
#
# (%20 caps concurrent tasks - adjust to your cluster's fair-use policy.) To smoke-test a
# single task first, instead run:
#     sbatch --array=1-1 $PIPELINE_DIR/workflow/step04_interproscan/interproscan_array_job.slurm.sh
#
set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── ADAPT: the paths below ───
PROTEIN_DIR="$STEP04_OUT/proteins"
MANIFEST="$STEP04_OUT/interpro_manifest.tsv"
OUT_DIR="$STEP04_OUT/interproscan_results"
# ─────────────────────────────────────────────────────────────────────────────

[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST (run build_interpro_manifest.slurm.sh first)" >&2; exit 1; }

### ADAPT ### — put InterProScan (and the Java it needs) on the PATH.
#              `set +u` around the module calls: lmod's init script references
#              variables that may be unset, which would trip `set -u`. Restored
#              to `set -u` right after.
set +u
module load adoptopenjdk/11.0.12+7
module load interproscan/5.66-98.0
set -u
### END ADAPT ###

command -v interproscan.sh >/dev/null || { echo "ERROR: interproscan.sh not on PATH" >&2; exit 1; }

mkdir -p "$OUT_DIR"

# Look up this task's (genome, chunk_num, start_seq, end_seq) from the manifest by task_id
# (column 1) matching SLURM_ARRAY_TASK_ID - NOT by line number, since the manifest has a
# header row and this stays correct even if rows are ever reordered.
row=$(awk -F'\t' -v tid="$SLURM_ARRAY_TASK_ID" '$1 == tid {print; found=1} END {if (!found) exit 1}' "$MANIFEST") \
  || { echo "ERROR: no manifest row for task_id $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

genome=$(echo "$row"    | cut -f2)
chunk_num=$(echo "$row" | cut -f3)
start_seq=$(echo "$row" | cut -f4)
end_seq=$(echo "$row"   | cut -f5)

echo "=== Task $SLURM_ARRAY_TASK_ID: $genome, chunk $chunk_num (seqs $start_seq-$end_seq) ==="

fasta="$PROTEIN_DIR/${genome}.faa"
genome_out_dir="$OUT_DIR/$genome"
mkdir -p "$genome_out_dir"

chunk_file="$genome_out_dir/chunk_${chunk_num}.fasta"

awk -v start="$start_seq" -v end="$end_seq" '
BEGIN {seq_count=0}
/^>/ {
    seq_count++
    print_this = (seq_count >= start && seq_count <= end)
}
print_this {print}
' "$fasta" > "$chunk_file"

if [[ ! -s "$chunk_file" ]]; then
  echo "No sequences in this chunk, exiting"
  rm -f "$chunk_file"
  exit 0
fi

interproscan.sh \
    -i "$chunk_file" \
    -o "$genome_out_dir/chunk_${chunk_num}_interproscan.tsv" \
    -f tsv \
    -cpu "$SLURM_CPUS_PER_TASK" \
    -goterms \
    -iprlookup \
    -pa \
    -dp

rm -f "$chunk_file"
echo "Finished $genome chunk $chunk_num"
