#!/usr/bin/env bash
#
# Create the folders a project needs before its first sbatch. Safe to re-run: anything that
# already exists is left alone.
#
#     bash setup_project.sh [project folder]     # default: this checkout (the pipeline is cloned
#                                                #          into the project, as on main)
#
# Makes genomes/ and gff_originals/ (the inputs), plus every log folder a job's
# #SBATCH --output/--error line writes into (logs/slurm/, expanded_annotation/logs/, ...).
# SLURM opens the log before the job starts, so a job can't create its own log folder; without it
# the job fails and leaves no log.

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
project="${1:-$PIPELINE_DIR}"

# The folder part of every "#SBATCH --output=<dir>/<file>" / "--error=" line in the pipeline's jobs.
mapfile -t log_dirs < <(
    grep -rhoE '^#SBATCH --(output|error)=[^ ]+/' --include='*.sh' --include='*.slurm' \
        --exclude-dir=results --exclude-dir=.git --exclude-dir=.venv "$PIPELINE_DIR" |
        sed -E 's/^#SBATCH --(output|error)=//; s#/$##' | sort -u
)

mkdir -p "$project"
for dir in genomes gff_originals "${log_dirs[@]}"; do
    mkdir -p "$project/$dir"
done

echo "project ready: $(cd "$project" && pwd)"
echo "  genomes/        put the assemblies here (<name>.fna)"
echo "  gff_originals/  put each genome's own annotation here (<name>.gff3)"
echo "  log folders:    ${log_dirs[*]}"
