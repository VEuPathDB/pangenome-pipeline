#!/bin/bash
#SBATCH --job-name=liftoff_all_vs_all
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=200:00:00
#
# All-against-all Liftoff: every genome's annotation is lifted onto every other genome's FASTA.
#
# Input annotation = the GFFs written by prepare_for_liftoff_job.sh ($STEP05_OUT/gff_for_liftoff/),
# with liftoff_feature_types.txt from the same folder. FASTAs = the project's genomes/ folder.
#
# Submit from your project folder (logs/ must exist there - SLURM does not create it). Either
# everything in one job (sequential):
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/liftoff_all_vs_all.sh
# or one job per source genome, in parallel (N = number of genomes in GENOMES):
#     N=$(grep -c . genome_order.txt)
#     sbatch --array=1-${N} $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/liftoff_all_vs_all.sh
# Logs: logs/liftoff_all_vs_all_<jobid>.out / .err - one pair per array task.
#
# Safe to resubmit: a source->target pair whose polished output already exists is skipped.
# Parallel-safe: every pair uses its own temporary folder ($STEP05_OUT/liftoff_output/tmp/<source>_onto_<target>,
# deleted when the pair finishes; a folder left behind means that pair failed - see the .err).

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── Edit these for your run ───
GENOMES="genome_order.txt"                    # file with one genome name per line
GFF_DIR="$STEP05_OUT/gff_for_liftoff"        # prepared GFFs: <genome>.gff3 (or .gff)
FEATURE_TYPES="$GFF_DIR/liftoff_feature_types.txt"   # top-level types to transfer (written by prepare_for_liftoff.py).
                                              # Liftoff only transfers type "gene" unless a type is listed here (-f) -
                                              # exact match - so this is what keeps non-coding RNA genes from being copied.
FASTA_DIR="genomes"                           # <genome>.fna (or .fa / .fasta)
OUT_DIR="$STEP05_OUT/liftoff_output"         # results: ${OUT_DIR}/polished/<source>_onto_<target>_liftoff.gff3_polished
# ────────────────────────────────

set +u
module purge
module load liftoff/1.6.3
set -u

[[ -f "$GENOMES" ]]        || { echo "ERROR: genome list not found: $GENOMES (run from the project folder)" >&2; exit 1; }
[[ -f "$FEATURE_TYPES" ]]  || { echo "ERROR: $FEATURE_TYPES not found - run prepare_for_liftoff_job.sh first" >&2; exit 1; }

# Find "<dir>/<genome><ext>" for the first extension that exists.
find_file() {
    local dir="$1" genome="$2"; shift 2
    local ext
    for ext in "$@"; do
        [[ -f "$dir/${genome}${ext}" ]] && { echo "$dir/${genome}${ext}"; return 0; }
    done
    return 1
}

genomes=()
while IFS= read -r line; do
    line="${line%%#*}"; line="${line//[[:space:]]/}"
    [[ -n "$line" ]] && genomes+=("$line")
done < "$GENOMES"

# Check every input up front - a missing file should stop the run now, not 100 hours in.
missing=0
for g in "${genomes[@]}"; do
    find_file "$GFF_DIR" "$g" .gff3 .gff >/dev/null          || { echo "ERROR: no GFF for $g in $GFF_DIR" >&2; missing=1; }
    find_file "$FASTA_DIR" "$g" .fna .fa .fasta >/dev/null   || { echo "ERROR: no FASTA for $g in $FASTA_DIR" >&2; missing=1; }
done
[[ $missing -eq 0 ]] || exit 1
echo "Found ${#genomes[@]} genomes: ${genomes[*]}"

# As a SLURM array, task N handles only the Nth genome as the source.
sources=("${genomes[@]}")
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    idx=$((SLURM_ARRAY_TASK_ID - 1))
    [[ $idx -ge 0 && $idx -lt ${#genomes[@]} ]] || { echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} is outside 1-${#genomes[@]}" >&2; exit 1; }
    sources=("${genomes[$idx]}")
    echo "Array task ${SLURM_ARRAY_TASK_ID}: source genome ${sources[0]}"
fi

mkdir -p "$OUT_DIR/polished" "$OUT_DIR/tmp"

for source in "${sources[@]}"; do
    for target in "${genomes[@]}"; do
        [[ "$source" == "$target" ]] && continue

        final="$OUT_DIR/polished/${source}_onto_${target}_liftoff.gff3_polished"
        if [[ -s "$final" ]]; then
            echo "Skipping ${source} onto ${target} - already done"
            continue
        fi

        echo "Lifting ${source} annotation onto ${target} genome"

        src_gff=$(find_file "$GFF_DIR" "$source" .gff3 .gff)
        src_fa=$(find_file "$FASTA_DIR" "$source" .fna .fa .fasta)
        tgt_fa=$(find_file "$FASTA_DIR" "$target" .fna .fa .fasta)

        # Liftoff writes temporary files (some with fixed names, e.g. polish.sam) into -dir, which
        # defaults to ./intermediate_files - shared by every run in the same folder. Runs in parallel
        # (array tasks) would overwrite each other, so every pair gets its own.
        tmpdir="$OUT_DIR/tmp/${source}_onto_${target}"
        rm -rf "$tmpdir"

        liftoff -g "$src_gff" \
                -o "$OUT_DIR/${source}_onto_${target}_liftoff.gff3" \
                -u "$OUT_DIR/${source}_onto_${target}_unmapped.txt" \
                -p 16 \
                -dir "$tmpdir" \
                -polish \
                -f "$FEATURE_TYPES" \
                -copies -sc 0.9 \
                -flank 0.1 \
                -mismatch 1 -gap_open 1 -gap_extend 1 \
                "$tgt_fa" \
                "$src_fa"

        mv "$OUT_DIR/${source}_onto_${target}_liftoff.gff3_polished" "$OUT_DIR/polished/"
        rm -rf "$tmpdir"
    done
done

echo "Liftoff completed for: ${sources[*]}. Polished outputs are in $OUT_DIR/polished/"
