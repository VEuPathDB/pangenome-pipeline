#!/bin/bash
#SBATCH --job-name=update_interpro_gff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Writes each genome's own step-04 InterProScan results into its own GFF (see
# update_interpro_in_gff.py for exactly how) - replacing any InterPro/PFAM evidence the
# original annotation already had, so every genome ends up on the same InterProScan version.
# One genome at a time, no cross-genome logic - that version of the problem is step 05's own,
# later, after Liftoff.
#
# Submit from your project folder:
#     sbatch $PIPELINE_DIR/workflow/step04_interproscan/update_interpro_in_gff.slurm.sh
#
set -eo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# ─── ADAPT: the paths below and the #SBATCH lines ───
GFF_INDIR="$STEP02_OUT/gff_cleaned"           # step 02's TE-kill-list output
LOOKUP_DIR="$STEP04_OUT/interpro_lookup"      # from merge_interpro_job.slurm.sh
ORDER="genome_order.txt"                       # step-01 canonical genome list
GFF_OUTDIR="$STEP04_OUT/gff_with_interpro"    # new GFFs land here - GFF_INDIR is untouched
# ────────────────────────────────────────────────────────────────────────────

### ADAPT ### — put python3 on the PATH:
activate_env "$ENV_BIOINFO"  # env names / micromamba: pipeline_config.sh
### END ADAPT ###

set -u

GFF_DIR="$GFF_INDIR"
LOOKUP_PATH="$LOOKUP_DIR"
ORDER_FILE="$ORDER"
OUT_DIR="$GFF_OUTDIR"

[[ -d "$GFF_DIR" ]]    || { echo "ERROR: GFF folder not found: $GFF_DIR" >&2; exit 1; }
[[ -d "$LOOKUP_PATH" ]] || { echo "ERROR: lookup folder not found: $LOOKUP_PATH (run merge_interpro_job.slurm.sh first)" >&2; exit 1; }
[[ -f "$ORDER_FILE" ]] || { echo "ERROR: genome list not found: $ORDER_FILE" >&2; exit 1; }
mkdir -p "$OUT_DIR"

gff_for () {  # basename -> path
    for e in gff3 gff; do
        [[ -f "$GFF_DIR/$1.$e" ]] && { echo "$GFF_DIR/$1.$e"; return 0; }
    done
    return 1
}

mapfile -t LABELS < <(grep -vE '^[[:space:]]*(#|$)' "$ORDER_FILE" | awk '{print $1}')

n=0
for label in "${LABELS[@]}"; do
    gff=$(gff_for "$label") || { echo "skip $label (no GFF in $GFF_DIR)"; continue; }
    lookup="$LOOKUP_PATH/${label}_interpro.tsv"
    if [[ ! -f "$lookup" ]]; then
        echo "skip $label (no $lookup - run merge_interpro_job.slurm.sh for it first)"
        continue
    fi
    ext="${gff##*.}"
    echo "=== $label ==="
    python3 -m pangenome_pipeline.step04_interproscan.update_interpro_in_gff --gff "$gff" --lookup "$lookup" --out "$OUT_DIR/${label}.${ext}"
    n=$((n + 1))
done
(( n > 0 )) || { echo "ERROR: no genome had both a GFF and a lookup table - nothing done" >&2; exit 1; }

echo
echo "Done — $n genome(s) updated, output in $OUT_DIR/"
