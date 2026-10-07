#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=drop_contigs
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=drop_contigs_%j.out
#SBATCH --error=drop_contigs_%j.err
#
# Remove named contigs (e.g. organelle contigs from the organelle screen) from
# the assemblies in genomes/ and the matching features from gff_originals/,
# driven by a  <genome> <contig>  list file. See 01_Genome_Selection_QC/README.md
# section 9.
#
# This is a seconds-long text job - it does not need a cluster; `bash` is fine.
# Submit / run from the project folder:
#     sbatch 01_Genome_Selection_QC/scripts/drop_contigs_slurm.sh
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#  1. LIST below - the path to your  <genome> <contig>  file.
#  2. The  #SBATCH  lines above (partition / resources).
#  3. The  ### ADAPT ###  line - it only needs a  python3  on the PATH. Use your
#     cluster's python module / `source activate <env>`, or delete it if python3
#     is already there.
# ═══════════════════════════════════════════════════════════════════════════════

# --- the list of contigs to drop -----------------------------------------------
LIST="01_Genome_Selection_QC/drop_contigs.txt"
# -----------------------------------------------------------------------------

set -eo pipefail

### ADAPT ### — put python3 on the PATH:
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###

set -u

PROJECT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
SCRIPT="$PROJECT_DIR/01_Genome_Selection_QC/scripts/drop_contigs.py"
LIST_PATH="$PROJECT_DIR/$LIST"

[[ -f "$SCRIPT" ]]    || { echo "ERROR: $SCRIPT not found" >&2; exit 1; }
[[ -f "$LIST_PATH" ]] || { echo "ERROR: list file $LIST_PATH not found - set LIST at the top" >&2; exit 1; }

cd "$PROJECT_DIR"

# 1. dry run - show exactly what will be removed
echo "=== DRY RUN ==="
python3 "$SCRIPT" "$LIST_PATH"

# 2. apply
echo
echo "=== APPLYING ==="
python3 "$SCRIPT" "$LIST_PATH" --apply

echo
echo "Done. genomes/ and gff_originals/ rewritten. Re-run the section 5.1 GFF"
echo "cross-check and the organelle screen to confirm the set is now clean."
