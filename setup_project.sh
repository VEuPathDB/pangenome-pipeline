#!/usr/bin/env bash
#
# Create a species project folder with the folders the pipeline expects. Safe to re-run: anything
# that already exists is left alone.
#
#     bash setup_project.sh /path/to/my_species_project
#
# Then put the assemblies in genomes/ and their GFFs in gff_originals/, and submit jobs from inside
# the project folder (cd /path/to/my_species_project && sbatch $PIPELINE_DIR/workflow/...).

set -euo pipefail

[[ $# -eq 1 ]] || { echo "usage: bash $0 /path/to/my_species_project" >&2; exit 1; }

mkdir -p "$1"/{genomes,gff_originals,logs,output}
echo "project ready: $(cd "$1" && pwd)"
echo "  genomes/        put the assemblies here (<name>.fna)"
echo "  gff_originals/  put each genome's own annotation here (<name>.gff3)"
echo "  logs/           SLURM logs"
echo "  output/         every step's results"
