#!/bin/bash
#SBATCH --job-name=te_gene_overlap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=te_gene_overlap_%j.out
#SBATCH --error=te_gene_overlap_%j.err
#
# Full TE/gene-overlap review in one job, all output in one folder (OUTDIR):
#   1. te_gene_overlap.py       - per-gene TE overlap -> gene_te_overlap.tsv (gene_id,
#                                 gene_description, coordinates, TE_overlap_pct,
#                                 TE_divergence_pct - the RepeatMasker %-divergence-from-
#                                 consensus of the TE fragment that explains most of the
#                                 gene's overlap, etc. - all the per-gene info lives in this
#                                 one table), plots/, protein_coding_threshold_summary.tsv,
#                                 pseudogene_threshold_summary.tsv
#   2-3. extract_flagged_genes.py + extract_cds_sequences.py, run three times for three
#        non-overlapping bands (GFF records + spliced CDS FASTA for each):
#          - "review"           MIN_PCT <= overlap <  MAX_PCT             (inspect band -
#                                mandatory check before eliminating any of these)
#          - "eliminated"                overlap >= MAX_PCT               (all of them)
#          - "eliminated_review"         overlap >= MAX_PCT  AND  divergence >= REVIEW_MIN_DIV
#                                (eliminated but on a divergent/fragmented TE call - possible
#                                domesticated gene or spurious call, worth a manual look
#                                before trusting the elimination)
#
# Submit from your project folder:
#     sbatch 02_TE_annotation_Clean_Up/scripts/te_gene_overlap.slurm.sh
#
# ─── ADAPT: the paths below, the ### ADAPT ### block, and the #SBATCH lines ───
TE_INDIR="."                          # folder with the RM_<genome_id>/ dirs (each with one .out)
GFF_INDIR="gff_originals"             # one GFF3 per genome
GENOMES="genomes"                     # step-01 assemblies folder (for the CDS extraction)
OUTDIR="te_annotation_integration"    # everything (tsv/plots/flagged GFFs/CDS fasta) goes here
GENOME_MAP=""                         # optional: genome_id<TAB>gff_path TSV (step 1 only)
MIN_PCT=20                            # "review" band lower bound (the "inspect" band)
MAX_PCT=80                            # "review" band upper bound == "eliminated" band lower bound
REVIEW_MIN_DIV=30                     # eliminated genes with TE_divergence_pct >= this get a
                                       # second, "eliminated_review" band flagged for manual check
EXTRA_GENES=""                        # optional: genome_label<TAB>gene_id TSV, no header
                                       # (added to the "review" band's CDS extraction only)
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script PY 02_TE_annotation_Clean_Up/scripts/te_gene_overlap.py
SCRIPTS="$(dirname "$PY")"

### ADAPT ### — put a python3 with pandas + matplotlib on the PATH
#              (the 'synteny' env from step 01 once pandas is added, or
#               `micromamba install -c conda-forge pandas matplotlib`):
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
python3 -c "import pandas, matplotlib" 2>/dev/null \
    || { echo "ERROR: python3 is missing pandas and/or matplotlib" >&2; exit 1; }
compgen -G "$TE_INDIR"/RM_*/ >/dev/null || compgen -G "$TE_INDIR"/R_*/ >/dev/null \
    || compgen -G "$TE_INDIR"/*.out >/dev/null \
    || { echo "ERROR: no RM_<genome_id>/ folders or *.out files in $TE_INDIR" >&2; exit 1; }
[[ -d "$GFF_INDIR" ]] || { echo "ERROR: GFF folder not found: $GFF_INDIR" >&2; exit 1; }
[[ -d "$GENOMES" ]]   || { echo "ERROR: genomes folder not found: $GENOMES" >&2; exit 1; }

echo "=== 1/4: te_gene_overlap.py ==="
OVERLAP_ARGS=(--te-indir "$TE_INDIR" --gff-indir "$GFF_INDIR" --outdir "$OUTDIR")
[[ -n "${GENOME_MAP:-}" ]] && OVERLAP_ARGS+=(--genome-map "$GENOME_MAP")
python3 "$SCRIPTS/te_gene_overlap.py" "${OVERLAP_ARGS[@]}"
GENE_TABLE="$OUTDIR/gene_te_overlap.tsv"

echo
echo "=== 2/4: review band (${MIN_PCT}% <= overlap < ${MAX_PCT}%) - mandatory check ==="
python3 "$SCRIPTS/extract_flagged_genes.py" \
    --gene-table "$GENE_TABLE" --gff-indir "$GFF_INDIR" --outdir "$OUTDIR" \
    --min-pct "$MIN_PCT" --max-pct "$MAX_PCT" --out-prefix flagged_genes
CDS_ARGS=(--gene-table "$GENE_TABLE" --gff-indir "$GFF_INDIR" --genomes-dir "$GENOMES" \
          --outdir "$OUTDIR" --min-pct "$MIN_PCT" --max-pct "$MAX_PCT" \
          --out-name flagged_genes_cds.fna)
[[ -n "${EXTRA_GENES:-}" ]] && CDS_ARGS+=(--extra-genes "$EXTRA_GENES")
python3 "$SCRIPTS/extract_cds_sequences.py" "${CDS_ARGS[@]}"

echo
echo "=== 3/4: eliminated band (overlap >= ${MAX_PCT}%) - all of them ==="
python3 "$SCRIPTS/extract_flagged_genes.py" \
    --gene-table "$GENE_TABLE" --gff-indir "$GFF_INDIR" --outdir "$OUTDIR" \
    --min-pct "$MAX_PCT" --max-pct 100.01 --out-prefix eliminated_genes
python3 "$SCRIPTS/extract_cds_sequences.py" \
    --gene-table "$GENE_TABLE" --gff-indir "$GFF_INDIR" --genomes-dir "$GENOMES" \
    --outdir "$OUTDIR" --min-pct "$MAX_PCT" --max-pct 100.01 \
    --out-name eliminated_genes_cds.fna

echo
echo "=== 4/4: eliminated_review band (overlap >= ${MAX_PCT}%, divergence >= ${REVIEW_MIN_DIV}%) ==="
python3 "$SCRIPTS/extract_flagged_genes.py" \
    --gene-table "$GENE_TABLE" --gff-indir "$GFF_INDIR" --outdir "$OUTDIR" \
    --min-pct "$MAX_PCT" --max-pct 100.01 --min-div "$REVIEW_MIN_DIV" \
    --out-prefix eliminated_genes_review
python3 "$SCRIPTS/extract_cds_sequences.py" \
    --gene-table "$GENE_TABLE" --gff-indir "$GFF_INDIR" --genomes-dir "$GENOMES" \
    --outdir "$OUTDIR" --min-pct "$MAX_PCT" --max-pct 100.01 --min-div "$REVIEW_MIN_DIV" \
    --out-name eliminated_genes_review_cds.fna

echo
echo "Done — everything in $OUTDIR/:"
echo "  gene_te_overlap.tsv, plots/ (every gene, full table)"
echo "  flagged_genes.tsv / <genome_id>.flagged_genes.gff / flagged_genes_cds.fna"
echo "    (review band, ${MIN_PCT}-${MAX_PCT}%: mandatory check before eliminating)"
echo "  eliminated_genes.tsv / <genome_id>.eliminated_genes.gff / eliminated_genes_cds.fna"
echo "    (>= ${MAX_PCT}%: everything auto-eliminated)"
echo "  eliminated_genes_review.tsv / <genome_id>.eliminated_genes_review.gff /"
echo "  eliminated_genes_review_cds.fna"
echo "    (>= ${MAX_PCT}% but divergence >= ${REVIEW_MIN_DIV}%: double-check before trusting)"
