#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=check_species
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=check_species_%j.out
#SBATCH --error=check_species_%j.err
#
# Confirm every genome is the same species as the reference (line 1 of
# genome_order.txt), by whole-genome ANI with skani. Catches mislabelled public
# assemblies - especially cryptic sibling species - and gross contamination.
#
# Submit from the project folder:
#     sbatch 01_Genome_Selection_QC/scripts/check_species.slurm.sh
#
# Output:  species_check/species_check.tsv
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#  1. The  #SBATCH  lines above (partition / resources).
#  2. The TOOL ENVIRONMENT block (### ADAPT ###) — it only has to put  skani  on
#     the PATH. Replace with  module load skani /  source activate … / your paths.
#  3. THRESHOLD below, if the default 95% ANI is not right for your genus.
# ═══════════════════════════════════════════════════════════════════════════════

THRESHOLD=95        # ANI (%) at/above which a genome is called the same species

set -eo pipefail

### ADAPT ### — put skani on the PATH:
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###

set -u

command -v skani >/dev/null || {
    echo "ERROR: 'skani' is not on the PATH. Install it into your env, e.g.:" >&2
    echo "    micromamba install -n synteny -c bioconda -c conda-forge skani" >&2
    exit 1
}

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
GENOME_DIR="$PROJECT_DIR/genomes"
ORDER="$PROJECT_DIR/genome_order.txt"
OUT_DIR="$PROJECT_DIR/species_check"
mkdir -p "$OUT_DIR"

fasta_for () {
    for e in fna fasta fa; do
        [[ -f "$GENOME_DIR/$1.$e" ]] && { echo "$GENOME_DIR/$1.$e"; return 0; }
    done
    return 1
}

mapfile -t LABELS < <(grep -vE '^[[:space:]]*(#|$)' "$ORDER" | awk '{print $1}')
[[ ${#LABELS[@]} -ge 2 ]] || { echo "ERROR: need >= 2 genomes in $ORDER" >&2; exit 1; }

REF="${LABELS[0]}"
REF_FA=$(fasta_for "$REF") || { echo "ERROR: reference $REF not in $GENOME_DIR" >&2; exit 1; }

TSV="$OUT_DIR/species_check.tsv"
printf 'genome\tANI_to_%s\taligned_fraction\tverdict\n' "$REF" > "$TSV"
printf '%s\t%s\t%s\t%s\n' "$REF" "100.00" "1.00" "reference" >> "$TSV"

# skani dist columns: Ref_file Query_file ANI Align_fraction_ref Align_fraction_query ...
for label in "${LABELS[@]:1}"; do
    Q_FA=$(fasta_for "$label") || { echo "skip $label (no FASTA)"; continue; }
    out=$(skani dist -t "${SLURM_CPUS_PER_TASK:-4}" -q "$Q_FA" -r "$REF_FA") || out=""
    row=$(awk 'NR==2' <<< "$out")
    ani=$(awk '{print $3+0}' <<< "$row")
    af=$(awk '{printf "%.2f", ($5+0)/100}' <<< "$row")
    if [[ -z "$row" ]]; then           # skani reported nothing -> too divergent to align
        ani="0.00"; af="0.00"
    fi
    verdict=$(awk -v a="$ani" -v t="$THRESHOLD" 'BEGIN{print (a+0>=t)?"same_species":"CHECK"}')
    printf '%s\t%s\t%s\t%s\n' "$label" "$ani" "$af" "$verdict" >> "$TSV"
    echo "$label  ANI=$ani  af=$af  -> $verdict"
done

echo
echo "Genomes to review:"
awk -F'\t' '$4=="CHECK"{print "  "$1"  (ANI "$2")"}' "$TSV"
echo "Full table: $TSV"
