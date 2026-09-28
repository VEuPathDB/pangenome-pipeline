#!/bin/bash
#SBATCH --job-name=te_aid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --output=te_aid_%j.out
#SBATCH --error=te_aid_%j.err
#
# Run TE-Aid (https://github.com/clemgoub/TE-Aid) on every candidate TE consensus
# in a FASTA, to visually validate them one by one (coverage / structure / ORF /
# self dot-plot panel per candidate). TE-Aid takes ONE consensus per run, so this
# script splits the candidate FASTA and loops, pairing each candidate with the
# genome it was recovered from (TE-Aid's -g needs that genome).
#
# Submit from your project folder (the one that holds genomes/):
#     sbatch 02_TE_annotation_Clean_Up/scripts/run_te_aid.slurm.sh
#
# ─── ADAPT: the variables below, the ### ADAPT ### block, and the #SBATCH lines ──
#
# Candidate FASTA to inspect. Exactly one INPUT_FASTA / OUTDIR pair must be active
# (uncommented) at a time. To use a different input, comment the pair below and
# uncomment one of the others. Each pair writes its own OUTDIR, so you can run the
# script once per file without overwriting anything.
#
# Paths are relative to the folder you submit from (the project root — the one with
# genomes/). The candidate-pipeline outputs live under 02_TE_annotation_Clean_Up/,
# hence the prefix; adjust it if yours are elsewhere.
#
#   DEFAULT — the filtered set (only candidates with a RepBase / TE-peptide hit, Step 4 output):
INPUT_FASTA="02_TE_annotation_Clean_Up/final_TE_candidates_after_MCHelper/final_potential_new_TEs.fa"
OUTDIR="te_aid_final_potential_new_TEs"
#
#   the unfiltered non-redundant set — also shows the candidates the RepBase/peptide
#   filter dropped (e.g. solo-LTR-only families). Comment the pair above, uncomment this one:
# INPUT_FASTA="02_TE_annotation_Clean_Up/final_TE_candidates_after_MCHelper/all_MCHelper_curated_sequences_NR.cdhit80.fa"
# OUTDIR="te_aid_all_MCHelper_NR_cdhit80"
#
#   Complementary TE Screen pipeline output, if you ran it. Comment the DEFAULT pair, uncomment this one:
# INPUT_FASTA="02_TE_annotation_Clean_Up/TE_gap_screen/final_gap_screen_candidates.cdhit80.fa"
# OUTDIR="te_aid_gap_screen_candidates"

GENOMES="genomes"                 # <genome_id><suffix> assemblies (step-01 folder), or a full path
GENOME_SUFFIX=".fna"              # removed from the assembly file name to get <genome_id>
TEAID="/home/aludwig/hc-storage/bin/TE-Aid/TE-Aid"   # path to the TE-Aid executable
# OUTDIR is set together with INPUT_FASTA above — one sub-folder per candidate is created under it
# ───────────────────────────────────────────────────────────────────────────────

set -euo pipefail

### ADAPT ### — put TE-Aid's dependencies on the PATH: R (with its TE-Aid packages),
#              NCBI BLAST+ (makeblastdb / blastn / tblastn) and hmmer.
#              `set +u` around the module/micromamba calls: lmod's init and the
#              micromamba hook reference variables that may be unset, which would
#              trip `set -u`. Restored to `set -u` right after.
set +u
module load ncbi-blast/2.12.0 || true
module load R/4.5.1 || true
source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
activate_env "$ENV_TE_ANNOT"  # env names / micromamba: pipeline_config.sh
set -u
### END ADAPT ###

[[ -s "$INPUT_FASTA" ]] || { echo "ERROR: candidate FASTA not found or empty: $INPUT_FASTA (run from the project folder)" >&2; exit 1; }
[[ -d "$GENOMES" ]]     || { echo "ERROR: genomes folder not found: $GENOMES" >&2; exit 1; }
[[ -x "$TEAID" ]]       || { echo "ERROR: TE-Aid not found or not executable: $TEAID" >&2; exit 1; }
command -v blastn >/dev/null || { echo "ERROR: blastn not on PATH (adapt the env block)" >&2; exit 1; }
command -v Rscript >/dev/null || { echo "ERROR: Rscript not on PATH (adapt the env block)" >&2; exit 1; }

PROJECT_DIR="$(pwd -P)"
mkdir -p "$OUTDIR/single_fastas"

echo "Candidate FASTA : $INPUT_FASTA"
echo "Genomes folder  : $GENOMES  (suffix '$GENOME_SUFFIX')"
echo "TE-Aid          : $TEAID"
echo "Output folder   : $OUTDIR"
echo

# ── split the multi-FASTA: one file per candidate, header sanitised for a filename ──
echo "=== Splitting $INPUT_FASTA ==="
awk -v outdir="$OUTDIR/single_fastas" '
    /^>/ {
        if (out) close(out)
        header = substr($0, 2)
        fname  = header
        gsub(/[^A-Za-z0-9_.-]/, "_", fname)
        out = outdir "/" fname ".fa"
        print ">" header > out
        next
    }
    { print $0 >> out }
' "$INPUT_FASTA"

# ── run TE-Aid per candidate, against the genome it came from ──
echo "=== Running TE-Aid ==="
summary="$OUTDIR/te_aid_run_summary.tsv"
printf "candidate_fasta\tgenome_id\tgenome_fasta\tstatus\n" > "$summary"

for fa in "$OUTDIR/single_fastas"/*.fa; do
    [[ -s "$fa" ]] || continue
    base="$(basename "$fa" .fa)"

    # Candidate headers from the pipeline are "<genome_id>|<candidate_id>"; the awk
    # step turned '|' into '_'. Take the leading token, then recover a full
    # GCA_/GCF_ accession if that is what the genome_id is.
    genome_id="${base%%_*}"
    acc="$(printf '%s\n' "$base" | grep -oE '^GC[AF]_[0-9]+\.[0-9]+' || true)"
    [[ -n "$acc" ]] && genome_id="$acc"

    genome="$GENOMES/${genome_id}${GENOME_SUFFIX}"
    if [[ ! -s "$genome" ]]; then
        echo "WARNING: genome not found for $base -> $genome" >&2
        printf "%s\t%s\t%s\tmissing_genome\n" "$fa" "$genome_id" "$genome" >> "$summary"
        continue
    fi

    dest="$OUTDIR/$base"
    mkdir -p "$dest"
    echo "  $base   vs   $genome_id"
    "$TEAID" -q "$PROJECT_DIR/$fa" -g "$PROJECT_DIR/$genome" -o "$PROJECT_DIR/$dest"
    printf "%s\t%s\t%s\tdone\n" "$fa" "$genome_id" "$genome" >> "$summary"
done

echo
echo "=== Done ==="
echo "One TE-Aid folder per candidate under $OUTDIR/ ; run table: $summary"
