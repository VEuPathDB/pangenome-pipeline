#!/bin/bash
#SBATCH --job-name=remove_genes_from_gff
#SBATCH --partition=low
#SBATCH --output=remove_genes_from_gff_%j.out
#SBATCH --error=remove_genes_from_gff_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#
# Removes the genes listed in a manual kill list (and everything nested under them) from each
# genome's GFF - see remove_genes_from_gff.py's own header for exactly what it does.
#
# Submit from your project folder:
#     sbatch 05_Basket_annotation_integration/scripts/remove_genes_from_gff_job.sh

set -euo pipefail

# ─── Edit these for your run ───
GENOMES="genome_order.txt"                         # comma-separated genome names, or a file with one per line
GFF_DIR="gff_with_interpro"                        # GFFs to remove genes from (never edited in place)
KILL_LIST="kill_list_manual.tsv"                   # TSV, no header: genome<TAB>gene_id
OUT_DIR="gff_curated"                              # where the curated GFFs go (every genome, kill-listed or not)
# ────────────────────────────────

set +u
module purge
### ADAPT ### — put micromamba on the PATH, then activate the env (standard library only):
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate bioinfo
### END ADAPT ###
set -u

# Locate the .py. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there - fall back to the standard layout
# relative to the directory you submitted from.
PY=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/remove_genes_from_gff.py" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Basket_annotation_integration/scripts/remove_genes_from_gff.py" \
    "05_Basket_annotation_integration/scripts/remove_genes_from_gff.py"
do
    [[ -f "$cand" ]] && { PY="$cand"; break; }
done
[[ -n "$PY" ]] || { echo "ERROR: cannot find remove_genes_from_gff.py (run from the project folder)" >&2; exit 1; }

[[ -f "$KILL_LIST" ]] || { echo "ERROR: kill list not found: $KILL_LIST" >&2; exit 1; }

python3 "$PY" \
  --genomes "$GENOMES" --gff-dir "$GFF_DIR" \
  --kill-list "$KILL_LIST" --out-dir "$OUT_DIR"

echo "Done. Curated GFFs are in $OUT_DIR"
