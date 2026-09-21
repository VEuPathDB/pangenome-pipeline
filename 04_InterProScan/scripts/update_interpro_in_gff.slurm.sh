#!/bin/bash
#SBATCH --job-name=update_interpro_gff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=update_interpro_gff_%j.out
#SBATCH --error=update_interpro_gff_%j.err
#
# Writes each genome's own step-04 InterProScan results into its own GFF (see
# update_interpro_in_gff.py for exactly how) - replacing any InterPro/PFAM evidence the
# original annotation already had, so every genome ends up on the same InterProScan version.
# One genome at a time, no cross-genome logic - that version of the problem is step 05's own,
# later, after Liftoff.
#
# Submit from your project folder:
#     sbatch 04_InterProScan/scripts/update_interpro_in_gff.slurm.sh
#
# ─── ADAPT: the paths below and the #SBATCH lines ───
GFF_INDIR="gff_cleaned"                        # step 02's TE-kill-list output
LOOKUP_DIR="interpro/interpro_lookup"          # from merge_interpro_job.slurm.sh
ORDER="genome_order.txt"                       # step-01 canonical genome list
GFF_OUTDIR="gff_with_interpro"                 # new GFFs land here - GFF_INDIR is untouched
# ────────────────────────────────────────────────────────────────────────────

set -eo pipefail

### ADAPT ### — put python3 on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate bioinfo
### END ADAPT ###

set -u

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
GFF_DIR="$PROJECT_DIR/$GFF_INDIR"
LOOKUP_PATH="$PROJECT_DIR/$LOOKUP_DIR"
ORDER_FILE="$PROJECT_DIR/$ORDER"
OUT_DIR="$PROJECT_DIR/$GFF_OUTDIR"
PY="$PROJECT_DIR/04_InterProScan/scripts/update_interpro_in_gff.py"

[[ -f "$PY" ]]         || { echo "ERROR: $PY not found" >&2; exit 1; }
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
    python3 "$PY" --gff "$gff" --lookup "$lookup" --out "$OUT_DIR/${label}.${ext}"
    n=$((n + 1))
done
(( n > 0 )) || { echo "ERROR: no genome had both a GFF and a lookup table - nothing done" >&2; exit 1; }

echo
echo "Done — $n genome(s) updated, output in $OUT_DIR/"
