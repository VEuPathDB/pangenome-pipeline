#!/bin/bash
#SBATCH --job-name=protein_evidence
#SBATCH --partition=low
#SBATCH --output=reference_check/logs/%x_%j.out
#SBATCH --error=reference_check/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00

set -euo pipefail

# Step 07, reference gene check (optional): protein evidence for the genes the pipeline added to the reference.
# protein_evidence.py searches each new protein with NCBI blastp against nr over the internet (blastp -remote), leaving
# out the reference's own species, and writes OUT_DIR/protein_evidence.tsv (best hit outside the species, or no hit).
# Only OUT_DIR is written. NCBI is slow: hundreds of proteins can take many hours. If the job stops (time limit,
# network), submit it again: the batches already done are kept and only the rest is searched.
# The compute node needs internet access; the job checks it first and stops if there is none.
# Then rebuild the page with the evidence (reference_check_job.sh with PROTEIN_EVIDENCE=OUT_DIR/protein_evidence.tsv).
# Submit from your project folder; the log folder must exist BEFORE submitting:
#     mkdir -p reference_check/logs
#     sbatch 07_Post_GET_PANGENES_adjusting/scripts/protein_evidence_job.sh

# ===== Edit these =====  (each can also be set on the command line: sbatch --export=ALL,NAME=value ...)
GENOME_ORDER="${GENOME_ORDER-genome_order.txt}"                                   # the reference is its first line
GENOME="${GENOME-}"                                                               # or name the reference here
STEP05_DIR="${STEP05_DIR-basket_annotation}"                                      # step 05 results
GFF_DIR="${GFF_DIR-$STEP05_DIR/gff_prefer_panel_model_final}"                     # step 05 final GFFs
GENOME_DIR="${GENOME_DIR-genomes}"                                                # <genome>.fna
EXCLUDE_TAXID="${EXCLUDE_TAXID-}"                                                 # the species left out (746128 A. fumigatus, 5811 T. gondii)
OUT_DIR="${OUT_DIR-reference_check/protein_evidence}"
# ======================

set +u
module purge
### ADAPT ### - python3 and blastp (BLAST+), as in step 05
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

SD=""
for cand in "$(dirname "${BASH_SOURCE[0]}")" \
            "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts" \
            "07_Post_GET_PANGENES_adjusting/scripts"; do
    [[ -f "$cand/protein_evidence.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find protein_evidence.py (run from the project folder)" >&2; exit 1; }

[[ -n "$EXCLUDE_TAXID" ]] || { echo "ERROR: set EXCLUDE_TAXID (746128 A. fumigatus, 5811 T. gondii)" >&2; exit 1; }
command -v blastp >/dev/null || { echo "ERROR: blastp not found (see ADAPT)" >&2; exit 1; }
[[ -n "$GENOME" ]] || GENOME=$(grep -vE '^\s*(#|$)' "$GENOME_ORDER" | head -1)
GFF="$GFF_DIR/${GENOME}_basket_annotation.gff"
FA=$(ls "$GENOME_DIR/$GENOME".f*a* 2>/dev/null | head -1 || true)
[[ -f "$GFF" ]] || { echo "ERROR: no $GFF" >&2; exit 1; }
[[ -n "$FA" ]] || { echo "ERROR: no $GENOME_DIR/$GENOME.fna" >&2; exit 1; }
if ! curl -s -m 30 -o /dev/null "https://blast.ncbi.nlm.nih.gov/Blast.cgi"; then
    echo "ERROR: no internet access from this node (blast.ncbi.nlm.nih.gov); remote BLAST is not possible here" >&2
    exit 1
fi
blastp -version | head -1
echo "reference: $GENOME ($GFF, $FA); species left out: txid$EXCLUDE_TAXID; output: $OUT_DIR"

python3 "$SD/protein_evidence.py" --gff "$GFF" --fasta "$FA" --exclude-taxid "$EXCLUDE_TAXID" --out-dir "$OUT_DIR"
echo "Done: $OUT_DIR/protein_evidence.tsv"
