#!/bin/bash
#SBATCH --job-name=RepeatMasker
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=RepeatMasker_%j.out
#SBATCH --error=RepeatMasker_%j.err
#
# Annotate TEs in every assembly with RepeatMasker, using the curated TE library.
# One output folder per genome: RM_<genome_id>/  (genome_id = the .fna file name
# without the extension; an NCBI accession is trimmed to GCA_/GCF_).
#
# Submit from your project folder (the one that holds genomes/):
#     sbatch 02_TE_annotation_Clean_Up/scripts/run_repeatmasker.slurm.sh
#
# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
GENOMES="genomes"                     # <genome_id>.fna assemblies (step-01 folder), or a full path
LIB="curated_TE_library.fa"           # curated TE library FASTA (in this folder), or a full path
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
### ADAPT ### — put RepeatMasker on the PATH (must be configured with a Dfam
#              library; see the README if `famdb.py … FamDB data directory not found`):
activate_env "$ENV_REPEATMASKER"  # env names / backend: pipeline_config.sh
### END ADAPT ###

command -v RepeatMasker >/dev/null || { echo "ERROR: RepeatMasker not on PATH" >&2; exit 1; }
[[ -d "$GENOMES" ]] || { echo "ERROR: genomes folder not found: $GENOMES (run from the project folder)" >&2; exit 1; }
[[ -s "$LIB" ]]     || { echo "ERROR: TE library not found: $LIB" >&2; exit 1; }

CPUS="${SLURM_CPUS_PER_TASK:-24}"
PA=$(( CPUS / 4 )); (( PA > 0 )) || PA=1        # rmblast uses 4 cores per -pa unit

for fna in "$GENOMES"/*.fna; do
    [[ -e "$fna" ]] || { echo "ERROR: no *.fna in $GENOMES" >&2; exit 1; }

    name=$(basename "$fna" .fna)
    acc=$(echo "$name" | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' || true)
    [[ -n "$acc" ]] && name="$acc"

    outdir="RM_${name}"
    mkdir -p "$outdir"

    echo "=========================================================="
    echo "Processing: $(basename "$fna")  ->  $outdir/   (-pa $PA)"
    echo "=========================================================="

    RepeatMasker \
        -pa "$PA" \
        -lib "$LIB" \
        -dir "$outdir" \
        -no_is \
        -gff \
        -xsmall \
        -html \
        "$fna"

    echo "Finished $name at $(date)"
done

echo "All RepeatMasker jobs completed successfully!"
