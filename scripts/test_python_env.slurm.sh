#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=test_python_env
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=test_python_env_%j.out
#SBATCH --error=test_python_env_%j.err
#
# Quick check that a compute node can run the pipeline's Python: the uv venv, its packages, the
# shared helpers in lib/pangenome_common, and (if CONTAINER_DIR is set) a container tool next to it.
#
# Submit from the project folder (the pipeline checkout), like every other job:
#     sbatch scripts/test_python_env.slurm.sh

set -eo pipefail

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # activate_env, pipeline_config.sh
activate_env "$ENV_SYNTENY"

set -u

echo "node:          $(hostname)"
echo "PIPELINE_DIR:  $PIPELINE_DIR"
echo "PIPELINE_VENV: $PIPELINE_VENV"
echo "python3:       $(command -v python3)"

[[ "$(command -v python3)" == "$PIPELINE_VENV/bin/python3" ]] \
    || echo "WARNING: python3 is not the pipeline venv's - run 'uv sync --frozen --no-dev' in $PIPELINE_DIR" >&2

PYTHONPATH="$PIPELINE_DIR/lib" python3 - <<'EOF'
import sys
import matplotlib, numpy, pandas, scipy, seaborn
import pangenome_common.gff

print("python       ", sys.version.split()[0])
for mod in (matplotlib, numpy, pandas, scipy, seaborn):
    print(f"{mod.__name__:<13}", mod.__version__)
print("lib           ", pangenome_common.gff.__file__)

# Plotting must work without a display on a compute node.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1])
fig.savefig("/dev/null", format="png")
print("matplotlib    plot OK")
EOF

python3 "$PIPELINE_DIR/01_Genome_Selection_QC/scripts/plot_dotplot.py" --help >/dev/null
echo "script --help OK"

if command -v minimap2 >/dev/null 2>&1; then
    echo "minimap2      $(minimap2 --version)  ($(command -v minimap2))"
fi

echo "ALL OK"
