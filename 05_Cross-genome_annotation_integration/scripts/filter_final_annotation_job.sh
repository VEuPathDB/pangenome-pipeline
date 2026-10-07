#!/bin/bash
#SBATCH --job-name=filter_final_annotation
#SBATCH --partition=low
#SBATCH --output=filtered_annotation/logs/%x_%j.out
#SBATCH --error=filtered_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00

set -euo pipefail

# Submit from your project folder. The log folder must exist BEFORE submitting (SLURM does not
# create it - without it the job leaves no log):
#     mkdir -p filtered_annotation/logs && sbatch 05_Cross-genome_annotation_integration/scripts/filter_final_annotation_job.sh

# Edit these for your run.
STRAINS="genome_order.txt"                 # comma-separated names or a file with one per line, first = the reference (same as ORDER)
FINAL_DIR="expanded_annotation"            # build_final_annotation_basic.sh's OUTPUT_DIR - the raw, unfiltered <genome>_final_annotation.gff files
OUT_DIR="filtered_annotation"              # where to write the filtered GFFs + reports
MIN_RESCUED_COVERAGE="0.5"                 # every rescued_<genome> gene below this Liftoff coverage is dropped outright
MIN_OVERLAP_BP="20"                        # CDS-to-CDS overlaps under this many bp are ignored - boundary-rounding noise
MIN_OVERLAP_GENOMES="2"                    # a gene present in >2 genomes needs a qualifying overlap in at least this many of them to count as a real, recurring conflict
EXTRA_COPY_REQUIRE_VALID_ORF="true"        # true/false - drop an extra-copy gene that has no valid ORF
REQUIRE_DOMAIN_ONLY="false"                # true = only a protein domain counts as evidence; false = a specific description also counts
ID_PREFIX="Afu"                            # renamed gene IDs look like <ID_PREFIX><N>g<00001, ...>

set +u
module purge
### ADAPT ### — put python3 on the PATH (standard library only, nothing else needed):
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###
set -u

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script PY 05_Cross-genome_annotation_integration/scripts/filter_final_annotation.py
mkdir -p "$OUT_DIR"

EXTRA_ARGS=()
[[ "$EXTRA_COPY_REQUIRE_VALID_ORF" == "false" ]] && EXTRA_ARGS+=(--no-extra-copy-require-valid-orf)
[[ "$REQUIRE_DOMAIN_ONLY" == "true" ]] && EXTRA_ARGS+=(--require-domain-only)

python3 "$PY" \
  --strains "$STRAINS" --final-dir "$FINAL_DIR" \
  --out-dir "$OUT_DIR" \
  --min-rescued-coverage "$MIN_RESCUED_COVERAGE" \
  --min-overlap-bp "$MIN_OVERLAP_BP" \
  --min-overlap-genomes "$MIN_OVERLAP_GENOMES" \
  --id-prefix "$ID_PREFIX" \
  ${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}

echo "All genomes filtered. Outputs are in $OUT_DIR"
