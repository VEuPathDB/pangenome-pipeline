#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=contig_redundancy
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=contig_redundancy_%j.out
#SBATCH --error=contig_redundancy_%j.err
#
# Report-only contig-redundancy screen (step 01). For every genome in
# genome_order.txt it:
#   1. indexes the FASTA          -> contig_redundancy/<label>.fai
#   2. self-aligns it             -> contig_redundancy/<label>.self.paf
#                                    (minimap2 -xasm5 -DP <asm> <asm>, the call purge_dups uses)
# then contig_redundancy.py combines that with the existing
# synteny/<label>_vs_<reference>.paf (reference = first genome in genome_order.txt, as in
# synteny_vs_ref.slurm.sh) and gff_originals/<label>.gff into
#   contig_redundancy/REPORT.md   +   contig_redundancy/<label>.contigs.tsv
#
# Nothing is dropped. Review the report, then edit
# 01_Genome_Selection_QC/drop_contigs.txt by hand and run drop_contigs_slurm.sh.
#
# Submit from the project folder:
#     sbatch 01_Genome_Selection_QC/scripts/contig_redundancy.slurm.sh
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#  1. The  #SBATCH  lines above.
#  2. The TOOL ENVIRONMENT block (### ADAPT ###) — it only has to put
#     minimap2, samtools  and  python3  on the PATH. As written it activates the
#     same micromamba env "synteny" used by synteny_vs_ref.slurm.sh.
#  3. The TUNING block — in particular CHROM_MIN_MB. The defaults suit an
#     A. fumigatus-sized genome; check them against YOUR species.
# ═══════════════════════════════════════════════════════════════════════════════

set -eo pipefail

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
### ADAPT ### — put minimap2 / samtools / python3 on the PATH:
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###

# ── TUNING ───────────────────────────────────────────────────────────────────
# CHROM_MIN_MB: contigs >= this size are treated as chromosome-scale ("large")
#   and are never flagged as redundant. Set it a little BELOW your species'
#   SMALLEST real chromosome / replicon, so a whole small chromosome is never
#   mistaken for a redundant fragment — and so a genuinely small extra contig is
#   still examined. 1 Mb fits most filamentous fungi (smallest A. fumigatus
#   chromosome ~1.8 Mb). Lower it for compact protist genomes; raise it if the
#   smallest chromosome is large. Look up your reference's chromosome sizes first.
CHROM_MIN_MB="1.0"
DUP_FRAC="90"       # % of a contig duplicated in a longer contig -> "small-redundant"
DUP_MIN_ID="0.95"   # min identity of the self-alignment blocks that count
# ─────────────────────────────────────────────────────────────────────────────

set -u

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
GENOME_DIR="$PROJECT_DIR/genomes"
ORDER="$PROJECT_DIR/genome_order.txt"
OUT_DIR="$PROJECT_DIR/contig_redundancy"
PY="$PROJECT_DIR/01_Genome_Selection_QC/scripts/contig_redundancy.py"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

[[ -f "$ORDER" ]] || { echo "ERROR: $ORDER not found" >&2; exit 1; }
[[ -f "$PY" ]]    || { echo "ERROR: $PY not found" >&2; exit 1; }
[[ -d "$GENOME_DIR" ]] || { echo "ERROR: $GENOME_DIR not found" >&2; exit 1; }

mkdir -p "$OUT_DIR"

fasta_for () {  # label -> path (.fna / .fasta / .fa)
    for e in fna fasta fa; do
        [[ -f "$GENOME_DIR/$1.$e" ]] && { echo "$GENOME_DIR/$1.$e"; return 0; }
    done
    return 1
}

mapfile -t LABELS < <(grep -vE '^[[:space:]]*(#|$)' "$ORDER" | awk '{print $1}')
[[ ${#LABELS[@]} -ge 1 ]] || { echo "ERROR: no genomes in $ORDER" >&2; exit 1; }

for label in "${LABELS[@]}"; do
    fa=$(fasta_for "$label") || { echo "skip $label (no FASTA in $GENOME_DIR)"; continue; }
    echo "=== $label ==="

    # 1. contig lengths
    samtools faidx "$fa"
    cp -f "${fa}.fai" "$OUT_DIR/${label}.fai"

    # 2. self-alignment (skip the self-diagonal with -D, keep all chains with -P)
    minimap2 -c -x asm5 -DP -t "$THREADS" "$fa" "$fa" > "$OUT_DIR/${label}.self.paf"
done

echo
echo "=== building report ==="
python3 "$PY" "$PROJECT_DIR" --ref "${LABELS[0]}" --md "contig_redundancy/REPORT.md" \
    --chrom-min-mb "$CHROM_MIN_MB" --dup-frac "$DUP_FRAC" --dup-min-id "$DUP_MIN_ID"

echo
echo "Done. See contig_redundancy/REPORT.md and contig_redundancy/<label>.contigs.tsv"
echo "Review, then add contigs to 01_Genome_Selection_QC/drop_contigs.txt by hand."
