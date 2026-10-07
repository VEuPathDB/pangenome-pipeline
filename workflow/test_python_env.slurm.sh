#!/bin/bash
#SBATCH --partition=low
#SBATCH --job-name=test_python_env
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# Quick check that a compute node can run the pipeline's Python: the uv venv, its packages, the
# pangenome_pipeline package, and (if CONTAINER_DIR is set) a container tool next to it.
#
# Submit from the project folder (logs/ must exist):
#     mkdir -p logs && sbatch $PIPELINE_DIR/workflow/test_python_env.slurm.sh

set -eo pipefail

_pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
[[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }
activate_env "$ENV_SYNTENY"

set -u

echo "node:          $(hostname)"
echo "PIPELINE_DIR:  $PIPELINE_DIR"
echo "PIPELINE_VENV: $PIPELINE_VENV"
echo "python3:       $(command -v python3)"

[[ "$(command -v python3)" == "$PIPELINE_VENV/bin/python3" ]] \
    || echo "WARNING: python3 is not the pipeline venv's - run 'uv sync --frozen --no-dev' in $PIPELINE_DIR" >&2

python3 - <<'EOF'
import sys
import matplotlib, numpy, pandas, scipy, seaborn
import pangenome_pipeline

print("python       ", sys.version.split()[0])
for mod in (matplotlib, numpy, pandas, scipy, seaborn):
    print(f"{mod.__name__:<13}", mod.__version__)
print("package       ", pangenome_pipeline.__path__[0])

# Plotting must work without a display on a compute node.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot([0, 1], [0, 1])
fig.savefig("/dev/null", format="png")
print("matplotlib    plot OK")
EOF

python3 -m pangenome_pipeline.step01_genome_selection_qc.plot_dotplot --help >/dev/null
echo "module --help OK"

if command -v minimap2 >/dev/null 2>&1; then
    echo "minimap2      $(minimap2 --version)  ($(command -v minimap2))"
fi

echo "ALL OK"
