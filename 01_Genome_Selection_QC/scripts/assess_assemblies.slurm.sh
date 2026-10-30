#!/bin/bash
#SBATCH --job-name=assess_assemblies
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=assess_assemblies_%j.out
#SBATCH --error=assess_assemblies_%j.err
#
# Quick assembly and annotation assessment of every genome: size, contigs, N50, GC, gene / mRNA /
# exon / CDS / intron counts and protein lengths (assess_assemblies.py), and a 7-panel figure
# comparing the genomes (plot_assemblies.py).
#
# Output, all in OUTDIR: assembly_metrics.tsv (one row per genome), protein_lengths.tsv (every
# protein), assembly_assessment.png.
#
# Submit from your project folder (the one that holds genomes/):
#     sbatch 01_Genome_Selection_QC/scripts/assess_assemblies.slurm.sh
#
# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
GENOMES="genomes"                       # step-01 assemblies folder (<name>.fna); NCBI's
                                        # <acc>_protein.faa / <acc>_data_summary.tsv may sit here too
GFF_INDIR="gff_originals"               # one GFF3 per genome, same basename as its FASTA
OUTDIR="assembly_assessment"            # where the tables/plot go
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# Locate the .py files. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there — fall back to the standard layout
# relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/01_Genome_Selection_QC/scripts" \
    "01_Genome_Selection_QC/scripts"
do
    [[ -f "$cand/assess_assemblies.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find assess_assemblies.py (run from the project folder)" >&2; exit 1; }

[[ -d "$GENOMES" ]] || { echo "ERROR: genomes folder not found: $GENOMES (run from the project folder)" >&2; exit 1; }
[[ -d "$GFF_INDIR" ]] || { echo "ERROR: GFF folder not found: $GFF_INDIR" >&2; exit 1; }

set +u
### ADAPT ### — put a python3 with pandas + matplotlib on the PATH:
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
python3 -c "import pandas, matplotlib" 2>/dev/null \
    || { echo "ERROR: python3 is missing pandas and/or matplotlib" >&2; exit 1; }

echo "=== 1/2: assess_assemblies.py ==="
python3 "$SD/assess_assemblies.py" --input-dir "$GENOMES" --gff-dir "$GFF_INDIR" --output-dir "$OUTDIR"

echo
echo "=== 2/2: plot_assemblies.py ==="
python3 "$SD/plot_assemblies.py" --metrics "$OUTDIR/assembly_metrics.tsv" \
    --proteins "$OUTDIR/protein_lengths.tsv" --output "$OUTDIR/assembly_assessment.png"

echo
echo "Done — $OUTDIR/assembly_metrics.tsv, $OUTDIR/protein_lengths.tsv, $OUTDIR/assembly_assessment.png"
