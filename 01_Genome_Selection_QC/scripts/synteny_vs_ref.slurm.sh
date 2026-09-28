#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=synteny_vs_ref
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=synteny_vs_ref_%j.out
#SBATCH --error=synteny_vs_ref_%j.err
#
# Whole-genome alignment (minimap2) of the reference against every other genome
# in genome_order.txt, plus a dot-plot per genome.
#
# Submit from the project folder:
#     sbatch 01_Genome_Selection_QC/scripts/synteny_vs_ref.slurm.sh
#
# Output:  synteny/<label>_vs_<ref>.paf   and   synteny/<label>_vs_<ref>.png
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING — two parts are specific to your cluster / setup:
#
#  1. The  #SBATCH  lines above — partition name, cpus, mem, time.
#
#  2. The TOOL ENVIRONMENT block below (### ADAPT ###). It only has to put
#     minimap2, samtools  and a  python3  with  matplotlib  on the PATH. As
#     written it activates a micromamba env called "synteny"; replace it with
#     whatever your site uses, e.g.  `module load minimap2 samtools`  /
#     `source activate synteny`  / your own micromamba paths.
# ═══════════════════════════════════════════════════════════════════════════════

set -eo pipefail

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
### ADAPT ### — put minimap2 / samtools / (python + matplotlib) on the PATH:
activate_env "$ENV_SYNTENY"  # env names / micromamba: pipeline_config.sh
### END ADAPT ###

set -u

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
GENOME_DIR="$PROJECT_DIR/genomes"
ORDER="$PROJECT_DIR/genome_order.txt"
OUT_DIR="$PROJECT_DIR/synteny"

PRESET="asm10"      # minimap2 -x : asm5 (<1% divergence) | asm10 | asm20 (more divergent)
MIN_ALN=2000        # drop alignments shorter than this from the plot

mkdir -p "$OUT_DIR"

mapfile -t LABELS < <(grep -vE '^[[:space:]]*(#|$)' "$ORDER" | awk '{print $1}')
[[ ${#LABELS[@]} -ge 2 ]] || { echo "ERROR: need >= 2 genomes in $ORDER" >&2; exit 1; }

fasta_for () {  # label -> path, .fna or .fasta
    for e in fna fasta fa; do
        [[ -f "$GENOME_DIR/$1.$e" ]] && { echo "$GENOME_DIR/$1.$e"; return 0; }
    done
    return 1
}

REF="${LABELS[0]}"
REF_FA=$(fasta_for "$REF") || { echo "ERROR: reference $REF not in $GENOME_DIR" >&2; exit 1; }
echo "Reference: $REF ($REF_FA)"

for label in "${LABELS[@]:1}"; do
    Q_FA=$(fasta_for "$label") || { echo "skip $label (no FASTA)"; continue; }
    paf="$OUT_DIR/${label}_vs_${REF}.paf"
    png="$OUT_DIR/${label}_vs_${REF}.png"
    echo "=== $label ==="
    minimap2 -x "$PRESET" -t "${SLURM_CPUS_PER_TASK:-8}" --secondary=no \
        "$REF_FA" "$Q_FA" > "$paf"
    python3 -m pangenome_pipeline.step01_genome_selection_qc.plot_dotplot "$paf" --ref-name "$REF" --query-name "$label" \
        --min-aln "$MIN_ALN" -o "$png"
done

echo "Done. See $OUT_DIR/"
