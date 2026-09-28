#!/bin/bash
#SBATCH --job-name=tag_final_annotation
#SBATCH --partition=low
#SBATCH --output=tagged_annotation/logs/%x_%j.out
#SBATCH --error=tagged_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00

set -euo pipefail

# Alternative to filter_final_annotation_job.sh: tags every gene with what it overlaps instead
# of dropping anything - see tag_final_annotation.py's docstring for why. Submit from your
# project folder. The log folder must exist BEFORE submitting (SLURM does not create it):
#     mkdir -p tagged_annotation/logs && sbatch 05_Cross-genome_annotation_integration/scripts/tag_final_annotation_job.sh

# Edit these for your run.
STRAINS="genome_order.txt"                 # comma-separated names or a file with one per line, first = the reference (same as ORDER)
FINAL_DIR="expanded_annotation"            # build_final_annotation_basic.sh's OUTPUT_DIR - the raw, unfiltered <genome>_final_annotation.gff files
OUT_DIR="tagged_annotation"                # where to write the tagged GFFs + intersect-pair reports
MIN_OVERLAP_BP="20"                        # CDS-to-CDS overlaps under this many bp are ignored - boundary-rounding noise
ID_PREFIX="Afu"                            # renamed gene IDs look like <ID_PREFIX><N>g<00001, ...> - edit this per species

set +u
module purge
### ADAPT ### — put python3 on the PATH (standard library only, nothing else needed):
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###
set -u

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script PY 05_Cross-genome_annotation_integration/scripts/tag_final_annotation.py
mkdir -p "$OUT_DIR"

python3 "$PY" \
  --strains "$STRAINS" --final-dir "$FINAL_DIR" \
  --out-dir "$OUT_DIR" \
  --min-overlap-bp "$MIN_OVERLAP_BP" \
  --id-prefix "$ID_PREFIX"

echo "All genomes tagged (nothing removed). Outputs are in $OUT_DIR"
