#!/bin/bash
# Run every pangenome plot for one pangene matrix. Runs locally (seconds, no SLURM needed).
#
#   plot_all.sh -i <pangene_matrix.tr.tab> -o <out_dir> -l <label>
#
# Example:
#   plot_all.sh -i output/step06_get_pangenes/for_get_pangenes_7b_pangenes/<subfolder>/pangene_matrix.tr.tab \
#       -o output/step08_pangenome_visualisation -l Toxo7b
#
# The category cuts can be changed with the environment variables CORE_MISSING, SOFTCORE_MISSING and CLOUD_MAX
# (default 0, 1 and 2 genomes); plot_all_job.sh sets them from its options.
# Every output file starts with the label, e.g. Toxo7b_pancore.png, Toxo7b_occupancy.svg.
# Needs python3 with numpy, scipy and matplotlib.

set -euo pipefail

MATRIX=""; OUT_DIR=""; LABEL=""
while getopts "i:o:l:h" opt; do
    case "$opt" in
        i) MATRIX="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        l) LABEL="$OPTARG" ;;
        *) sed -n '2,11p' "$0"; exit 1 ;;
    esac
done
[[ -n "$MATRIX" && -n "$OUT_DIR" && -n "$LABEL" ]] || { sed -n '2,11p' "$0"; exit 1; }
[[ -f "$MATRIX" ]] || { echo "ERROR: matrix not found: $MATRIX" >&2; exit 1; }

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# the pipeline's Python package, so this also works when run directly (not only from plot_all_job.sh)
export PYTHONPATH="$DIR/../../src${PYTHONPATH:+:$PYTHONPATH}"
mkdir -p "$OUT_DIR"

python3 -m pangenome_pipeline.step08_pangenome_visualisation.plot_pancore   -i "$MATRIX" -o "$OUT_DIR/${LABEL}_pancore"
python3 -m pangenome_pipeline.step08_pangenome_visualisation.plot_pancore   -i "$MATRIX" -o "$OUT_DIR/${LABEL}_pancore_separate" --layout separate
python3 -m pangenome_pipeline.step08_pangenome_visualisation.plot_occupancy -i "$MATRIX" -o "$OUT_DIR/${LABEL}_occupancy" \
    --core-missing "${CORE_MISSING:-0}" --softcore-missing "${SOFTCORE_MISSING:-1}" --cloud-max "${CLOUD_MAX:-2}"
python3 -m pangenome_pipeline.step08_pangenome_visualisation.plot_singletons -i "$MATRIX" -o "$OUT_DIR/${LABEL}_singletons"

# POCS.matrix.tab is written by get_pangenes next to pangene_matrix.tr.tab
POCS="$(dirname "$MATRIX")/POCS.matrix.tab"
if [[ -f "$POCS" ]]; then
    python3 -m pangenome_pipeline.step08_pangenome_visualisation.plot_pocs -i "$POCS" -o "$OUT_DIR/${LABEL}_POCS_heatmap"
else
    echo "WARNING: $POCS not found, skipping the POCS heatmap" >&2
fi

echo "Done. Figures for $LABEL are in $OUT_DIR"
