#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=screen_organelles
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=screen_organelles_%j.out
#SBATCH --error=screen_organelles_%j.err
#
# For each genome in genome_order.txt: is a given organelle present as its own
# contig, spread / DUPLICATED over several contigs, only a fragment inserted in a
# nuclear sequence, or absent?
#
# You tell it ONE genome that you know contains the organelle contig(s), and the
# contig name(s). It pulls those out, blastn's them against every genome, and
# writes one row per (genome, organelle, contig) that hits (so multi-copy cases
# are visible).
#
# Submit from the project folder:
#     sbatch 01_Genome_Selection_QC/scripts/screen_organelles.slurm.sh
#
# Output:  organelle_screen/organelle_screen.tsv  + a summary in the .out log
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#
#  1. The two variables under "WHICH ORGANELLE(S)" below.
#
#  2. The  #SBATCH  lines above (partition / resources).
#
#  3. The TOOL ENVIRONMENT block (### ADAPT ###) — it only has to put  blastn ,
#     samtools  and  python3  on the PATH. Replace with  module load … /
#     source activate … / your own micromamba paths.
# ═══════════════════════════════════════════════════════════════════════════════

# --- WHICH ORGANELLE(S) to look for ------------------------------------------
# SOURCE_GENOME     : basename (in genomes/) of a genome that you KNOW contains
#                     the organelle contig(s). Does not have to be the reference.
# ORGANELLE_CONTIGS : the contig name(s) in that genome, space-separated. You
#                     find these by looking at the assembly description and the
#                     contig names of that genome.
SOURCE_GENOME="${SOURCE_GENOME:?set SOURCE_GENOME (e.g. AF293) here or with sbatch --export=ALL,SOURCE_GENOME=...}"
ORGANELLE_CONTIGS="${ORGANELLE_CONTIGS:?set ORGANELLE_CONTIGS (e.g. mito_A_fumigatus_Af293) here or with sbatch --export}"
# --------------------------------------------------------------------------

set -eo pipefail

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
### ADAPT ### — put blastn / samtools / python3 on the PATH:
activate_env "$ENV_SYNTENY"  # env names / micromamba: pipeline_config.sh
### END ADAPT ###

set -u

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
GENOME_DIR="$PROJECT_DIR/genomes"
ORDER="$PROJECT_DIR/genome_order.txt"
OUT_DIR="$PROJECT_DIR/organelle_screen"
CLASSIFY="$PROJECT_DIR/01_Genome_Selection_QC/scripts/screen_organelles.py"

mkdir -p "$OUT_DIR/blast"

fasta_for () {  # basename -> path
    for e in fna fasta fa; do
        [[ -f "$GENOME_DIR/$1.$e" ]] && { echo "$GENOME_DIR/$1.$e"; return 0; }
    done
    return 1
}

# 1. pull the organelle contig(s) out of the source genome
SRC=$(fasta_for "$SOURCE_GENOME") \
    || { echo "ERROR: source genome '$SOURCE_GENOME' not in $GENOME_DIR" >&2; exit 1; }
ORGANELLES="$OUT_DIR/organelles.fna"
# shellcheck disable=SC2086
samtools faidx "$SRC" $ORGANELLE_CONTIGS > "$ORGANELLES"
echo "Organelle query: $(grep -c '^>' "$ORGANELLES") sequence(s) from $SOURCE_GENOME"

# 2. blastn every genome against it
mapfile -t LABELS < <(grep -vE '^[[:space:]]*(#|$)' "$ORDER" | awk '{print $1}')
FMT='6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'
for label in "${LABELS[@]}"; do
    fa=$(fasta_for "$label") || { echo "skip $label (no FASTA)"; continue; }
    echo "=== $label ==="
    blastn -query "$ORGANELLES" -subject "$fa" -dust no \
        -outfmt "$FMT" -evalue 1e-10 -max_target_seqs 100 \
        > "$OUT_DIR/blast/${label}.tsv" || true
done

# 3. classify
python3 "$CLASSIFY" --blast-dir "$OUT_DIR/blast" --organelles "$ORGANELLES" \
    --out "$OUT_DIR/organelle_screen.tsv"

echo "Done. See $OUT_DIR/organelle_screen.tsv"
