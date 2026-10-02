#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=check_redundant_genes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# For every `small-redundant` contig that carries genes (from the contig-redundancy
# screen), check whether those genes are genuinely already present elsewhere in the
# same assembly — the confirmation to run before adding the contig to drop_contigs.txt.
#
#   for each gene on the flagged contig:
#     samtools faidx  -> cut its genomic span
#     blastn          -> against the whole assembly
#     keep the best hit NOT on the gene's own contig
#     GFF check       -> does that hit land on an annotated gene?
#
# Output: $STEP01_OUT/contig_redundancy/<label>.redundant_genes.tsv  + a per-contig verdict in
# the job log (DROP-OK = every gene confirmed elsewhere; REVIEW = some not).
#
# Submit from the project folder:
#     sbatch $PIPELINE_DIR/workflow/step01_genome_selection_qc/check_redundant_genes.slurm.sh
#     # one genome only:
#     sbatch $PIPELINE_DIR/workflow/step01_genome_selection_qc/check_redundant_genes.slurm.sh GCA_040142975.1
#
# ═══════════════════════════════════════════════════════════════════════════════
#  ADAPT BEFORE RUNNING:
#  1. #SBATCH lines above.
#  2. TOOL ENVIRONMENT (### ADAPT ###) — needs samtools, blastn, python3 on PATH.
#     As written it activates the same micromamba env "synteny".
# ═══════════════════════════════════════════════════════════════════════════════

set -eo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }
### ADAPT ### — put samtools / blastn / python3 on the PATH:
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###

set -u


REDUNDANCY_DIR="$STEP01_OUT/contig_redundancy"   # contig_redundancy.slurm.sh's output, read and added to
[[ -d "$REDUNDANCY_DIR" ]] || {
    echo "ERROR: $REDUNDANCY_DIR not found - run contig_redundancy.slurm.sh first" >&2; exit 1; }

# any args passed to sbatch are genome labels; none -> auto-detect
python3 -m pangenome_pipeline.step01_genome_selection_qc.check_redundant_genes "$PROJECT_DIR" --redundancy-dir "$REDUNDANCY_DIR" "$@"

echo
echo "Done. Per genome in $REDUNDANCY_DIR/ :"
echo "  <label>.redundant_genes.tsv        per-gene table + descriptions + verdict"
echo "  <label>.redundant_genes.fa         the extracted gene spans (manual BLAST / TE-aid)"
echo "  <label>.redundant_genes.blast.tsv  the raw blastn hits"
echo "DROP-OK  -> every gene on that contig has a >=98% / >=90%-coverage copy on an"
echo "           annotated gene elsewhere: safe to add the contig to drop_contigs.txt."
echo "REVIEW   -> some genes not confirmed: inspect those rows before dropping."
