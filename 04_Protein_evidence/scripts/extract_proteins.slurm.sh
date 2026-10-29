#!/bin/bash
#SBATCH --job-name=extract_proteins
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=extract_proteins_%j.out
#SBATCH --error=extract_proteins_%j.err
#
# Extracts one protein FASTA per genome (translated CDS, spliced across exons)
# from GFF_INDIR + GENOMES using AGAT - the input InterProScan actually needs.
# Not part of the colleague's original InterProScan scripts (those assume the
# .faa files already exist); this is the missing first step.
#
# Submit from your project folder:
#     sbatch 04_Protein_evidence/scripts/extract_proteins.slurm.sh
#
# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
# GFF_INDIR - which GFFs to extract from:
#   - "gff_cleaned"   (default) - step 02's TE-kill-list output, for species that got a
#                       TE library and went through step 02.
#   - "gff_originals" - for a species with NO TE library (step 02 skipped). By this point
#                       gff_originals/ has already been edited in place by step 01's
#                       "Remove the unwanted contigs" (organelle + redundant-contig
#                       removal), so it's the right fallback, not the raw downloaded GFF.
GFF_INDIR="gff_cleaned"
GENOMES="genomes"                 # step-01 assemblies folder
OUTDIR="interpro/proteins"        # <genome_label>.faa land here
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

[[ -d "$GFF_INDIR" ]] || { echo "ERROR: GFF folder not found: $GFF_INDIR (run from the project folder)" >&2; exit 1; }
[[ -d "$GENOMES" ]]   || { echo "ERROR: genomes folder not found: $GENOMES" >&2; exit 1; }
mkdir -p "$OUTDIR"

### ADAPT ### — put AGAT (agat_sp_extract_sequences.pl) on the PATH:
export MAMBA_ROOT_PREFIX="/home/aludwig/hc-storage/micromamba"
export MAMBA_EXE="/home/aludwig/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate bioinfo
### END ADAPT ###

command -v agat_sp_extract_sequences.pl >/dev/null \
    || { echo "ERROR: agat_sp_extract_sequences.pl not on PATH" >&2; exit 1; }

n=0
for gff in "$GFF_INDIR"/*.gff "$GFF_INDIR"/*.gff3; do
    [[ -e "$gff" ]] || continue
    label=$(basename "$gff"); label="${label%.gff3}"; label="${label%.gff}"

    fasta=""
    for ext in fna fasta fa; do
        cand="$GENOMES/${label}.${ext}"
        [[ -f "$cand" ]] && { fasta="$cand"; break; }
    done
    if [[ -z "$fasta" ]]; then
        echo "WARNING: no genome FASTA found for $label in $GENOMES - skipping" >&2
        continue
    fi

    echo "=== $label ==="
    agat_sp_extract_sequences.pl \
        --gff "$gff" \
        --fasta "$fasta" \
        -p -t cds \
        -o "$OUTDIR/${label}.faa"

    # AGAT writes a trailing "*" for the stop codon on every translated protein.
    # InterProScan's own FASTA loader rejects "*" outright (crashes with a Java
    # exception in Protein.setSequenceAndMd5, not a clean error message) - strip
    # it from sequence lines (never from header lines) before InterProScan sees it.
    sed -i '/^>/!s/\*//g' "$OUTDIR/${label}.faa"

    n=$((n + 1))
done
(( n > 0 )) || { echo "ERROR: no *.gff/*.gff3 files in $GFF_INDIR" >&2; exit 1; }

echo
echo "Done — $n protein FASTA(s) in $OUTDIR/"
