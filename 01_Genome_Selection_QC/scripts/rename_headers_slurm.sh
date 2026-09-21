#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=rename_headers
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=rename_headers_%j.out
#SBATCH --error=rename_headers_%j.err
#
# Trim FASTA headers, in place, in genomes/ (see 01_Genome_Selection_QC/README.md).
# Submit from the project folder:
#     sbatch 01_Genome_Selection_QC/scripts/rename_headers_slurm.sh
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#  1. The  #SBATCH  lines above (partition / resources).
#  2. The  module load  line below (### ADAPT ###) — it only needs a  python3 .
#     Use your cluster's python module, or  source activate <env>, or delete the
#     line if python3 is already on the PATH.
# ═══════════════════════════════════════════════════════════════════════════════

set -eo pipefail

### ADAPT ### — put python3 on the PATH:
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###

set -u

# The folder sbatch was run from (SLURM runs the script from a spool dir, so we
# can't use the script's own path). Falls back to $PWD for a plain `bash` run.
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# --- edit these if your folders differ -------------------------------------
SCRIPT="$PROJECT_DIR/01_Genome_Selection_QC/scripts/rename_headers.py"
GENOME_DIR="$PROJECT_DIR/genomes"
GFF_DIR="$PROJECT_DIR/gff_originals"     # set to "" to skip the GFF ID cross-check
# --------------------------------------------------------------------------

[[ -f "$SCRIPT" ]] || { echo "ERROR: $SCRIPT not found" >&2; exit 1; }

shopt -s nullglob
inputs=( "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fa )
shopt -u nullglob
[[ ${#inputs[@]} -ge 1 ]] || { echo "ERROR: no FASTA files in $GENOME_DIR" >&2; exit 1; }

python3 "$SCRIPT" "${inputs[@]}" --in-place

if [[ -n "$GFF_DIR" ]]; then
    shopt -s nullglob
    fastas=( "$GENOME_DIR"/*.fna "$GENOME_DIR"/*.fasta "$GENOME_DIR"/*.fa )
    shopt -u nullglob
    python3 "$SCRIPT" "${fastas[@]}" --check-gff "$GFF_DIR"
fi

echo "Done. $GENOME_DIR/ headers are trimmed."
