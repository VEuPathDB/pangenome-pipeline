#!/bin/bash
#SBATCH --job-name=explore_multigene_pangenes
#SBATCH --partition=low
#SBATCH --output=multigene_pangenes/logs/%x_%j.out
#SBATCH --error=multigene_pangenes/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail

# Run this after get_pangenes has finished. Submit from your project folder. The log folder must
# exist BEFORE submitting (SLURM does not create it):
#     mkdir -p multigene_pangenes/logs && sbatch 07_Post_GET_PANGENES_adjusting/scripts/explore_multigene_pangenes_job.sh

# Edit these for your run.
PANGENES_OUT_DIR="for_get_pangenes_7b_pangenes"   # get_pangenes output folder (holds the *2neigh_algMmap_ subfolder)
TAGGED_DIR="tagged_annotation"                    # step 05's 7b OUT_DIR (tagged GFFs + <genome>_intersect_pairs.tsv)
OUT_DIR="multigene_pangenes"
LABEL="Toxo7b"                                    # prefix of every output file

set +u
module purge
### ADAPT ### — put python3 on the PATH (standard library only, nothing else needed):
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
module load python/3.11.9
### END ADAPT ###
set -u

# Locate the .py. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there - fall back to the standard layout
# relative to the directory you submitted from.
PY=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")/explore_multigene_pangenes.py" \
    "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts/explore_multigene_pangenes.py" \
    "07_Post_GET_PANGENES_adjusting/scripts/explore_multigene_pangenes.py"
do
    [[ -f "$cand" ]] && { PY="$cand"; break; }
done
[[ -n "$PY" ]] || { echo "ERROR: cannot find explore_multigene_pangenes.py (run from the project folder)" >&2; exit 1; }
mkdir -p "$OUT_DIR"

# get_pangenes writes pangene_matrix_genes.tr.tab inside <genome>_..._2neigh_algMmap_/ under its output folder.
shopt -s nullglob
MATRICES=("$PANGENES_OUT_DIR"/*2neigh_algMmap_/pangene_matrix_genes.tr.tab)
shopt -u nullglob
if [[ ${#MATRICES[@]} -ne 1 ]]; then
    echo "ERROR: expected exactly one $PANGENES_OUT_DIR/*2neigh_algMmap_/pangene_matrix_genes.tr.tab, found ${#MATRICES[@]}" >&2
    exit 1
fi

python3 "$PY" \
  --genes-matrix "${MATRICES[0]}" \
  --tagged-dir "$TAGGED_DIR" \
  -o "$OUT_DIR/${LABEL}_multigene"

echo "Multi-gene pangene exploration completed. Outputs are in $OUT_DIR"
