#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=drop_contigs
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Remove named contigs (e.g. organelle contigs from the organelle screen) from
# the assemblies in genomes/ and the matching features from gff_originals/,
# driven by a  <genome> <contig>  list file. See docs/step01_genome_selection_qc.md
# section 12.
#
# This is a seconds-long text job - it does not need a cluster; `bash` is fine.
# Submit / run from the project folder:
#     sbatch $PIPELINE_DIR/workflow/step01_genome_selection_qc/drop_contigs_slurm.sh
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#  1. LIST below - the path to your  <genome> <contig>  file.
#  2. The  #SBATCH  lines above (partition / resources).
#  3. The  ### ADAPT ###  line - it only needs a  python3  on the PATH. Use your
#     cluster's python module / `source activate <env>`, or delete it if python3
#     is already there.
# ═══════════════════════════════════════════════════════════════════════════════

set -eo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# --- the list of contigs to drop -----------------------------------------------
LIST="$STEP01_OUT/drop_contigs.txt"   # written by hand (see the doc), kept with the step-01 outputs
# -----------------------------------------------------------------------------

### ADAPT ### — put python3 on the PATH:
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###

set -u

LIST_PATH="$LIST"

[[ -f "$LIST_PATH" ]] || { echo "ERROR: list file $LIST_PATH not found - set LIST at the top" >&2; exit 1; }

cd "$PROJECT_DIR"

# 1. dry run - show exactly what will be removed
echo "=== DRY RUN ==="
python3 -m pangenome_pipeline.step01_genome_selection_qc.drop_contigs "$LIST_PATH"

# 2. apply
echo
echo "=== APPLYING ==="
python3 -m pangenome_pipeline.step01_genome_selection_qc.drop_contigs "$LIST_PATH" --apply

echo
echo "Done. genomes/ and gff_originals/ rewritten. Re-run the section 5.1 GFF"
echo "cross-check and the organelle screen to confirm the set is now clean."
