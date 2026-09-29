#!/bin/bash
#SBATCH --job-name=tag_final_annotation
#SBATCH --partition=low
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

# Alternative to filter_final_annotation_job.sh: tags every gene with what it overlaps instead
# of dropping anything - see tag_final_annotation.py's docstring for why. Submit from your
# project folder (logs/ must exist there - SLURM does not create it):
#     sbatch $PIPELINE_DIR/workflow/step05_cross_genome_annotation_integration/tag_final_annotation_job.sh

set -euo pipefail

# Find the pipeline (PIPELINE_DIR, else this checkout, else above the submit folder), then load
# lib/slurm_common.sh: PROJECT_DIR, OUTPUT_DIR / STEPNN_OUT, find_script, activate_env, pipeline_config.sh
_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }

# Edit these for your run.
STRAINS="genome_order.txt"                 # comma-separated names or a file with one per line, first = the reference (same as ORDER)
FINAL_DIR="$STEP05_OUT/expanded_annotation"   # build_final_annotation_basic.sh's FINAL_DIR - the raw, unfiltered <genome>_final_annotation.gff files
OUT_DIR="$STEP05_OUT/tagged_annotation"   # where to write the tagged GFFs + intersect-pair reports
MIN_OVERLAP_BP="20"                        # CDS-to-CDS overlaps under this many bp are ignored - boundary-rounding noise
ID_PREFIX="Afu"                            # renamed gene IDs look like <ID_PREFIX><N>g<00001, ...> - edit this per species

set +u
module purge
### ADAPT ### — put python3 on the PATH (standard library only, nothing else needed):
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###
set -u

mkdir -p "$OUT_DIR"

python3 -m pangenome_pipeline.step05_cross_genome_annotation_integration.tag_final_annotation \
  --strains "$STRAINS" --final-dir "$FINAL_DIR" \
  --out-dir "$OUT_DIR" \
  --min-overlap-bp "$MIN_OVERLAP_BP" \
  --id-prefix "$ID_PREFIX"

echo "All genomes tagged (nothing removed). Outputs are in $OUT_DIR"
