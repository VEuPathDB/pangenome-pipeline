#!/bin/bash
#SBATCH --job-name=psauron
#SBATCH --array=1-1       # one task per genome: always submit with --array=1-N on the sbatch command line (see below)
#SBATCH --partition=low
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=psauron/logs/psauron_%A_%a.out
#SBATCH --error=psauron/logs/psauron_%A_%a.err
#
# PSAURON on every protein of one genome per array task: the proteins written by extract_proteins.slurm.sh
# (the end-of-protein and internal "*" already removed there). Output: psauron/<genome>.csv (comment lines, then
# description,psauron_is_protein,in-frame_score; description = the protein ID, i.e. the mRNA ID). Step 05 reads this
# folder (PSAURON_DIR="psauron" in basket_annotation_job.sh) and maps the mRNA IDs to their genes.
#
# Submit from your project folder, after extract_proteins.slurm.sh. psauron/logs must exist BEFORE submitting:
#     mkdir -p psauron/logs
#     N=$(grep -v '^#' genome_order.txt | awk 'NF' | wc -l)
#     sbatch --array=1-${N} 04_Protein_evidence/scripts/psauron_array_job.slurm.sh
# A genome whose CSV is already complete (one row per protein) is skipped, so the same command resumes a run.

set -euo pipefail

# ===== Edit these =====
PROTEIN_DIR="interpro/proteins"     # <genome>.faa, from extract_proteins.slurm.sh
ORDER="genome_order.txt"            # one genome per line (task i = line i)
OUT_DIR="psauron"
### ADAPT ### — the PSAURON container (the group pipeline's, PSAURON 1.0.8; singularity must be on the PATH):
SIF="/mnt/hc-storage/groups/veupathdb_annotation/containers2/GAQA_psauron.sif"
### END ADAPT ###
# ======================

TASK="${SLURM_ARRAY_TASK_ID:-1}"
genome=$(grep -v '^#' "$ORDER" | awk 'NF{print $1}' | sed -n "${TASK}p")
[[ -n "$genome" ]] || { echo "ERROR: no genome on line $TASK of $ORDER" >&2; exit 1; }
faa="$PROTEIN_DIR/$genome.faa"
[[ -s "$faa" ]] || { echo "ERROR: protein FASTA not found or empty: $faa" >&2; exit 1; }
[[ -f "$SIF" ]] || { echo "ERROR: container not found: $SIF" >&2; exit 1; }
command -v singularity >/dev/null || { echo "ERROR: singularity is not on the PATH" >&2; exit 1; }
mkdir -p "$OUT_DIR"

out="$OUT_DIR/$genome.csv"
n=$(grep -c '^>' "$faa")
rows() { awk -F, 'h{n++} $1=="description"{h=1} END{print n+0}' "$1"; }
if [[ -s "$out" && "$(rows "$out")" -eq "$n" ]]; then
    echo "$genome: $out already has all $n proteins - skipped"
    exit 0
fi

echo "$genome: $n proteins | $faa -> $out | cpus: ${SLURM_CPUS_PER_TASK:-?}"
start=$(date +%s)
# real absolute paths, run from the real output folder, mounted: a container cannot follow a symlinked working
# directory (e.g. ~/hc-storage -> /mnt/hc-storage/...) and would not find relative paths
out_abs="$(realpath "$OUT_DIR")/$genome.csv"
mkdir -p "$OUT_DIR/running"
tmp_abs="$(realpath "$OUT_DIR/running")/$genome.csv"   # moved to OUT_DIR when complete (step 05 reads OUT_DIR/*.csv)
# headers cut to the protein ID (AGAT adds "gene=... seq_id=..." after it), so the scores carry the plain ID
in_abs="$(realpath "$OUT_DIR/running")/$genome.faa"
awk '/^>/{print $1; next} {print}' "$faa" > "$in_abs"
( cd "$(realpath "$OUT_DIR")" && \
  singularity exec -B "$(dirname "$out_abs")" "$SIF" psauron -i "$in_abs" -o "$tmp_abs" -p )
mv "$tmp_abs" "$out_abs"
rm -f "$in_abs"
end=$(date +%s)

got=$(rows "$out")
echo "$genome: done in $((end - start)) s; $got scores for $n proteins"
[[ "$got" -eq "$n" ]] || echo "WARNING: $genome: $((n - got)) proteins without a score (step 05 stops if one of them is a candidate)" >&2
