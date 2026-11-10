#!/bin/bash
#SBATCH --job-name=basket_annotation
#SBATCH --partition=low
#SBATCH --output=basket_annotation/logs/%x_%j.out
#SBATCH --error=basket_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00

set -euo pipefail

# Basket annotation (step 05, after the preparation and Liftoff): builds each genome's final GFF, one gene per panel gene, from the
# genomes' own annotations and the all-vs-all Liftoff transfers. Runs, in order:
#   0. prepare_pseudogenes.py      the reference's pseudogenes (no CDS) get their span as CDS, so the next steps use them
#   1. baskets_v2.py               baskets (same place, same strand), checked by place on the reference (then on the next
#                                  genomes), and their winners: the reference genes always win; a non-reference gene
#                                  wins with a domain, or with a valid ORF + >= MIN_AA aa + PSAURON >= the minimum
#                                  (it also finds the tandem arrays, tandem_arrays.tsv: their baskets keep each
#                                  genome's own copies)
#   2. decide_models.py            what each genome takes for each panel gene, in both modes (prefer_panel_model: the
#                                  valid model most similar to the panel gene, minimap2 + blastn; absent.tsv says why)
#   3. build_basket_annotation.py  the GFFs, once per mode in MODES (draft: before miniprot and the overlap rule)
#   4. repair_with_miniprot.py     genes with an invalid ORF: repaired with miniprot, or flagged (if RUN_MINIPROT=yes)
#   5. overlap_filter.py           the overlap rule (opposite strand), on the repaired GFFs -> the FINAL GFFs
#   6. check_basket_annotation.py  checks of the final GFFs
#   7. plot_basket_annotation.py   plots of the results, the modes side by side
#   8. gene_fate.py                every native gene: its final model, or the record of why it has none (stops with
#                                  an error when a gene has neither)
#   9. plot_gene_fate.py           from the gene fate: each genome's gene balance (own annotation -> final GFF) and the
#                                  genes eliminated by the quality rules, by their signs, and their PSAURON
#                                  histograms (plot_eliminated_psauron.py)
# Submit from your project folder. The log folder must exist BEFORE submitting:
#     mkdir -p basket_annotation/logs && sbatch 05_Basket_annotation_integration/scripts/basket_annotation_job.sh

# ===== Edit these for your species and run =====  (each can also be set on the command line: sbatch --export=ALL,NAME=value ...)
GENOMES="${GENOMES-genome_order.txt}"                           # one genome per line, the REFERENCE FIRST (its order gives the ID indexes)
NATIVE_DIR="${NATIVE_DIR-gff_for_liftoff}"                      # each genome's own annotation, as given to Liftoff (<genome>.gff3 / .gff)
TRANSFER_DIR="${TRANSFER_DIR-liftoff_output/polished}"          # <source>_onto_<target>_liftoff.gff3_polished (step 4)
DOMAIN_DIR="${DOMAIN_DIR-annotation_check}"                     # <genome>_gene_classification.tsv (has_domain of each gene): output of
                                                                #   analyze_input_gffs_job.sh (step 1), run on the GFFs with InterPro
GENOME_DIR="${GENOME_DIR-genomes}"                              # <genome><GENOME_SUFFIX>
GENOME_SUFFIX="${GENOME_SUFFIX-.fna}"
OLD_REFERENCE_GFF="${OLD_REFERENCE_GFF-}"                       # reference GFF BEFORE its manual curation, if there was one ("" = none)
CURATION_REMOVED_LIST="${CURATION_REMOVED_LIST-}"               # the curators' list of the genes they deleted (IDs of OLD_REFERENCE_GFF, one per line);
                                                                #   "" = every old gene whose ID is not in the reference GFF counts
ID_PREFIX="${ID_PREFIX-}"                                       # prefix of the new gene IDs, e.g. Tg -> Tg1g00001 (REQUIRED)
MODES="${MODES-keep_similar prefer_panel_model}"                # keep_similar: keep a genome's own model when similar; prefer_panel_model: the valid model most similar to the panel gene
# each genome's model (decide_models.py):
MIN_COV="${MIN_COV-0.7}"                                        # a valid model must cover >= this share of the panel gene's CDS
MIN_TCOV="${MIN_TCOV-0.7}"                                      # the panel gene's transfer without a valid ORF is kept (for miniprot) from this coverage
MIN_COV_RESCUE="${MIN_COV_RESCUE-0.5}"                          # rescue before absent: a valid model covering >= this share
MINIMAP2_MODULE="${MINIMAP2_MODULE-}"                           # module that puts minimap2 on the PATH ("" if it is already there)
# winners without a reference gene (non-reference genes without a domain):
MIN_AA="${MIN_AA-100}"                                          # shortest protein
PSAURON_DIR="${PSAURON_DIR-psauron}"                            # PSAURON scores of all proteins (step 04: psauron_array_job.slurm.sh); "" = not available
RUN_PSAURON="${RUN_PSAURON-no}"                                 # used only when PSAURON_DIR is "": yes = score the candidates here; no = no PSAURON rule
PSAURON_CMD="${PSAURON_CMD-singularity exec /mnt/hc-storage/groups/veupathdb_annotation/containers2/GAQA_psauron.sif psauron}"
PSAURON_MIN="${PSAURON_MIN-0.7}"                                # minimum PSAURON score when the reference does not have the sequence
PSAURON_MIN_IN_REF="${PSAURON_MIN_IN_REF-0.9}"                  # minimum when the gene's transfer onto the reference works
PROTECTED_GENES="${PROTECTED_GENES-}"                            # the shield: genome<TAB>gene<TAB>reason; a basket holding one never
                                                                #   loses for the quality rules and wins overlaps; "" = no shield
OVERLAP_MIN_STRAINS="${OVERLAP_MIN_STRAINS-2}"                  # overlap rule: removed when it loses in >= min(this, genomes where present)
RUN_MINIPROT="${RUN_MINIPROT-yes}"                              # yes / no
MINIPROT_MODULE="${MINIPROT_MODULE-miniprot/0.18}"              # module that puts miniprot on the PATH ("" if it is already there)
PLOT_LABEL="${PLOT_LABEL-panel}"                                # prefix of the plot files
# genes removed before step 05, for the gene balance plot (optional; several entries separated by "|", each FILE or
# FILE:REASON; see plot_gene_fate.py):
REMOVED_GENES="${REMOVED_GENES-}"                               # tables of genes (genome, gene), or {genome} GFFs, e.g. step 02's TE genes
REMOVED_CONTIGS="${REMOVED_CONTIGS-}"                           # tables of contigs (genome, contig) removed in step 01
ORIGINAL_GFF_DIR="${ORIGINAL_GFF_DIR-}"                         # the GFFs before the contigs were removed (the genes on REMOVED_CONTIGS)
OUT_DIR="${OUT_DIR-basket_annotation}"
THREADS="${THREADS-8}"                                          # match --cpus-per-task above
# ===============================================

[[ -n "$ID_PREFIX" ]] || { echo "ERROR: set ID_PREFIX at the top of this script" >&2; exit 1; }

set +u
module purge
### ADAPT ### — put blastn / makeblastdb / minimap2 / python3 on the PATH, and miniprot (PSAURON runs through PSAURON_CMD):
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
if [[ "$RUN_MINIPROT" == "yes" && -n "$MINIPROT_MODULE" ]]; then module load "$MINIPROT_MODULE"; fi
if [[ -n "$MINIMAP2_MODULE" ]]; then module load "$MINIMAP2_MODULE"; fi
### END ADAPT ###
set -u
for tool in python3 blastn blastp makeblastdb minimap2; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool is not on the PATH (see ### ADAPT ###)" >&2; exit 1; }
done

# Locate the scripts. Under `sbatch` the wrapper is copied to a spool dir, so
# ${BASH_SOURCE[0]} is useless there - fall back to the standard layout
# relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/05_Basket_annotation_integration/scripts" \
    "05_Basket_annotation_integration/scripts"
do
    [[ -f "$cand/baskets_v2.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find baskets_v2.py (run from the project folder)" >&2; exit 1; }
REF=$(grep -v '^#' "$GENOMES" | awk 'NF{print $1; exit}')
REF_GFF=""
for x in "$NATIVE_DIR/$REF.gff3" "$NATIVE_DIR/$REF.gff"; do [[ -f "$x" ]] && { REF_GFF="$x"; break; }; done
[[ -n "$REF_GFF" ]] || { echo "ERROR: no $NATIVE_DIR/$REF.gff3 or .gff" >&2; exit 1; }

PSAURON_ARGS=()
if [[ -n "$PSAURON_DIR" ]]; then
    ls "$PSAURON_DIR"/*.csv >/dev/null 2>&1 || { echo "ERROR: no *.csv in PSAURON_DIR=$PSAURON_DIR: run 04_Protein_evidence/scripts/psauron_array_job.slurm.sh first, or set PSAURON_DIR=\"\" and RUN_PSAURON=\"yes\"" >&2; exit 1; }
    PSAURON_ARGS=(--psauron-dir "$PSAURON_DIR")
elif [[ "$RUN_PSAURON" == "yes" ]]; then
    PSAURON_ARGS=(--run-psauron --psauron-cmd "$PSAURON_CMD")
fi
if [[ ${#PSAURON_ARGS[@]} -gt 0 ]]; then echo "PSAURON: ${PSAURON_ARGS[*]}"; else echo "PSAURON: off"; fi

mkdir -p "$OUT_DIR"
PROTECTED_ARGS=()
rm -f "$OUT_DIR/protected_baskets.tsv" "$OUT_DIR/protected_genes_saved.tsv"   # never reuse an older run's shield
if [[ -n "$PROTECTED_GENES" ]]; then
    [[ -f "$PROTECTED_GENES" ]] || { echo "ERROR: no PROTECTED_GENES=$PROTECTED_GENES" >&2; exit 1; }
    PROTECTED_ARGS=(--protected-genes "$PROTECTED_GENES")
    echo "shield: $PROTECTED_GENES"
else
    echo "shield: off"
fi
echo "== 0. the reference's pseudogenes made usable (CDS = their span; written back as gene + transcript + exon)"
python3 "$SD/prepare_pseudogenes.py" --in-dir "$NATIVE_DIR" --out-dir "$OUT_DIR/input/native" --reference "$REF"
python3 "$SD/prepare_pseudogenes.py" --in-dir "$TRANSFER_DIR" --out-dir "$OUT_DIR/input/transfers" --transfers --reference "$REF"
NATIVE_DIR="$OUT_DIR/input/native"
TRANSFER_DIR="$OUT_DIR/input/transfers"
REF_GFF="$NATIVE_DIR/$(basename "$REF_GFF")"
echo "== 1. baskets and winners"
python3 "$SD/baskets_v2.py" --native-dir "$NATIVE_DIR" --transfer-dir "$TRANSFER_DIR" --genomes "$GENOMES" \
  --genome-dir "$GENOME_DIR" --genome-suffix "$GENOME_SUFFIX" --domain-dir "$DOMAIN_DIR" \
  --min-aa "$MIN_AA" --psauron-min "$PSAURON_MIN" --psauron-min-in-ref "$PSAURON_MIN_IN_REF" \
  ${PSAURON_ARGS[@]+"${PSAURON_ARGS[@]}"} ${PROTECTED_ARGS[@]+"${PROTECTED_ARGS[@]}"} --threads "$THREADS" \
  --out-dir "$OUT_DIR"
echo "== 2. decisions"
python3 "$SD/decide_models.py" --baskets "$OUT_DIR/baskets_v2.tsv" --native-dir "$NATIVE_DIR" \
  --transfer-dir "$TRANSFER_DIR" --genome-dir "$GENOME_DIR" --genome-suffix "$GENOME_SUFFIX" --genomes "$GENOMES" \
  --min-cov "$MIN_COV" --min-tcov "$MIN_TCOV" --min-cov-rescue "$MIN_COV_RESCUE" --arrays "$OUT_DIR/tandem_arrays.tsv" \
  --threads "$THREADS" --out-dir "$OUT_DIR/decide"
for MODE in $MODES; do
    echo "== 3. draft GFFs, mode $MODE"
    python3 "$SD/build_basket_annotation.py" --baskets "$OUT_DIR/baskets_v2.tsv" \
      --decisions "$OUT_DIR/decide/decisions.tsv" --mode "$MODE" \
      --winners "$OUT_DIR/winners.tsv" --own-invalid "$OUT_DIR/decide/own_invalid_orf.tsv" \
      --array-models "$OUT_DIR/decide/array_models.tsv" --members "$OUT_DIR/members.tsv" \
      ${PROTECTED_GENES:+--protected-baskets "$OUT_DIR/protected_baskets.tsv"} \
      --native-dir "$NATIVE_DIR" --transfer-dir "$TRANSFER_DIR" --domain-dir "$DOMAIN_DIR" --genomes "$GENOMES" \
      ${OLD_REFERENCE_GFF:+--old-reference-gff "$OLD_REFERENCE_GFF"} \
      ${CURATION_REMOVED_LIST:+--curation-removed-list "$CURATION_REMOVED_LIST"} --id-prefix "$ID_PREFIX" \
      --out-dir "$OUT_DIR/gff_$MODE"
    REPAIRED="$OUT_DIR/gff_$MODE"
    if [[ "$RUN_MINIPROT" == "yes" ]]; then
        echo "== 4. miniprot, mode $MODE"
        python3 "$SD/repair_with_miniprot.py" --gff-dir "$OUT_DIR/gff_$MODE" --native-dir "$NATIVE_DIR" \
          --genome-dir "$GENOME_DIR" --genome-suffix "$GENOME_SUFFIX" --genomes "$GENOMES" --threads "$THREADS" \
          --out-dir "$OUT_DIR/gff_${MODE}_miniprot"
        REPAIRED="$OUT_DIR/gff_${MODE}_miniprot"
    fi
    FINAL="$OUT_DIR/gff_${MODE}_final"
    echo "== 5. overlap rule -> final GFFs, mode $MODE"
    python3 "$SD/overlap_filter.py" --gff-dir "$REPAIRED" --build-dir "$OUT_DIR/gff_$MODE" --genomes "$GENOMES" \
      --min-strains "$OVERLAP_MIN_STRAINS" --winners "$OUT_DIR/winners.tsv" \
      --out-dir "$FINAL"
    echo "== 6. checks, mode $MODE ($FINAL)"
    python3 "$SD/check_basket_annotation.py" --gff-dir "$FINAL" --genomes "$GENOMES" --reference-gff "$REF_GFF"
done

echo "== 7. plots"
python3 "$SD/plot_basket_annotation.py" --result-dir "$OUT_DIR" --genomes "$GENOMES" --native-dir "$NATIVE_DIR" \
  --modes "$MODES" --label "$PLOT_LABEL" --out-dir "$OUT_DIR/plots"

REMOVED_ARGS=()
IFS='|' read -r -a _rg <<< "$REMOVED_GENES"
for x in ${_rg[@]+"${_rg[@]}"}; do [[ -n "$x" ]] && REMOVED_ARGS+=(--removed-genes "$x"); done
IFS='|' read -r -a _rc <<< "$REMOVED_CONTIGS"
for x in ${_rc[@]+"${_rc[@]}"}; do [[ -n "$x" ]] && REMOVED_ARGS+=(--removed-contigs "$x"); done
[[ -n "$ORIGINAL_GFF_DIR" ]] && REMOVED_ARGS+=(--original-gff-dir "$ORIGINAL_GFF_DIR")
for MODE in $MODES; do
    echo "== 8. gene fate, mode $MODE"
    python3 "$SD/gene_fate.py" --native-dir "$NATIVE_DIR" --genomes "$GENOMES" --result-dir "$OUT_DIR" \
      --build-dir "$OUT_DIR/gff_$MODE" --final-dir "$OUT_DIR/gff_${MODE}_final" \
      --out "$OUT_DIR/gene_fate_$MODE.tsv"
    echo "== 9. gene balance and eliminated genes, mode $MODE"
    python3 "$SD/plot_gene_fate.py" --gene-fate "$OUT_DIR/gene_fate_$MODE.tsv" \
      --final-dir "$OUT_DIR/gff_${MODE}_final" --result-dir "$OUT_DIR" --reference "$REF" \
      --label "${PLOT_LABEL}_${MODE}" --min-aa "$MIN_AA" --psauron-min "$PSAURON_MIN" \
      --psauron-min-in-ref "$PSAURON_MIN_IN_REF" ${REMOVED_ARGS[@]+"${REMOVED_ARGS[@]}"} --out-dir "$OUT_DIR/plots"
    python3 "$SD/plot_eliminated_psauron.py" --eliminated "$OUT_DIR/plots/${PLOT_LABEL}_${MODE}_eliminated.tsv" \
      --label "${PLOT_LABEL}_${MODE}" --psauron-min "$PSAURON_MIN" --psauron-min-in-ref "$PSAURON_MIN_IN_REF" \
      --out-dir "$OUT_DIR/plots"
done
if grep -q $'\tNO RECORD\t' "$OUT_DIR"/gene_fate_*.tsv; then
    echo "ERROR: native genes without a final model and without a record (fate NO RECORD in $OUT_DIR/gene_fate_<mode>.tsv)" >&2
    exit 1
fi

echo "Done. Final GFFs: $OUT_DIR/gff_<mode>_final/<genome>_basket_annotation.gff; plots in $OUT_DIR/plots"
