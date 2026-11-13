#!/bin/bash
#SBATCH --job-name=build_explorer
#SBATCH --partition=low
#SBATCH --output=pangenome_explorer/logs/%x_%j.out
#SBATCH --error=pangenome_explorer/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail

# The interactive pangenome explorer (explorer.html + data/), and the same page packed in one file
# (explorer_standalone.html, opens with a double click, no server). Two uses:
#   preliminary (after the safe fix): the fixed matrix and step 05's GFFs, to look for genes to remove
#       mkdir -p pangenome_explorer/logs && sbatch --export=ALL,SET=preliminary 07_Post_GET_PANGENES_adjusting/scripts/build_explorer_job.sh
#   final (last part of step 07, default): the final pangenome
#       mkdir -p pangenome_explorer/logs && sbatch 07_Post_GET_PANGENES_adjusting/scripts/build_explorer_job.sh
# The output goes to $OUT_DIR/$SET. The log folder must exist BEFORE submitting.
# Optional layers (left out when their input is not there): each strain's own annotation under its genes (from
# step 05's gene fate and input GFFs: build_source_layer.py), the reference's VEuPathDB orthogroups (OG_TABLE), and
# the transposable elements of step 02's RepeatMasker runs (TE_DIR: build_te_layer.py).

# ===== Edit these =====
SET="${SET-final}"                                           # preliminary | final
TITLE="${TITLE-Pangenome Explorer}"
SPECIES="${SPECIES-}"                                        # e.g. Toxoplasma gondii (shown in italics)
GENOME_ORDER="${GENOME_ORDER-genome_order.txt}"              # the reference is its first line
SOURCE_DB="${SOURCE_DB-VEuPathDB}"                           # not shown any more (the page says "original"); kept for compatibility
SOURCE_URL="${SOURCE_URL-https://veupathdb.org/veupathdb/app/record/gene/}"   # its gene page: URL + gene ID
LINK_GENOMES="${LINK_GENOMES-}"                                        # strains annotated in that database (links only for them); "" = all
STEP05_DIR="${STEP05_DIR-basket_annotation}"                 # step 05 results: gene fate, input/native (own annotation)
GENOME_DIR="${GENOME_DIR-genomes}"                           # <genome>.fna: ORF of each own gene (optional)
OG_TABLE="${OG_TABLE-}"                                      # VEuPathDB "Genes by Ortholog Pattern" of the reference (optional)
TE_DIR="${TE_DIR-}"                                          # step 02 RepeatMasker results, one RM_<genome>/ each (optional)
FIXED_DIR="${FIXED_DIR-pangenome_fixed}"                     # preliminary: fix_pangenome_safe_job.sh output
STEP05_GFF_DIR="${STEP05_GFF_DIR-$STEP05_DIR/gff_prefer_panel_model_final}"  # preliminary: step 05 final GFFs
FINAL_DIR="${FINAL_DIR-pangenome_final}"                     # final: remove_panel_genes_job.sh output
WINNERS="${WINNERS-$STEP05_DIR/winners.tsv}"           # step 05: PSAURON of every panel gene (optional)
ENTRY_LIST="${ENTRY_LIST-$HOME/hc-storage/interpro_entry.list}"   # InterPro names (optional; IDs only when missing)
OUT_DIR="${OUT_DIR-pangenome_explorer}"
CORE_MISSING="${CORE_MISSING-0}"          # core = missing from at most this many genomes
SOFTCORE_MISSING="${SOFTCORE_MISSING-1}"  # soft-core = missing from at most this many genomes (and not core)
CLOUD_MAX="${CLOUD_MAX-2}"                # cloud = present in at most this many genomes; shell = the rest
# ======================

set +u
module purge 2>/dev/null || true
### ADAPT ### — put python3 on the PATH (standard library only):
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate synteny
### END ADAPT ###
set -u

# Locate the scripts. Under `sbatch` the wrapper is copied to a spool dir, so ${BASH_SOURCE[0]} is useless there -
# fall back to the standard layout relative to the directory you submitted from.
SD=""
for cand in \
    "$(dirname "${BASH_SOURCE[0]}")" \
    "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts" \
    "07_Post_GET_PANGENES_adjusting/scripts"
do
    [[ -f "$cand/build_explorer_data.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find build_explorer_data.py (run from the project folder)" >&2; exit 1; }

ABSENCES=""
if [[ "$SET" == preliminary ]]; then
    MATRIX="$FIXED_DIR/pangene_matrix_genes.tr.tab"
    GFF_DIR="$STEP05_GFF_DIR"
    [[ -f "$FIXED_DIR/absences_explained.tsv" ]] && ABSENCES="$FIXED_DIR/absences_explained.tsv"
    SUBTITLE="Preliminary pangenome (after the safe fix, before the removal of curated genes)"
elif [[ "$SET" == final ]]; then
    MATRIX="$FINAL_DIR/pangene_matrix_genes.tr.tab"
    GFF_DIR="$FINAL_DIR/gff"
    [[ -f "$FINAL_DIR/absences_explained.tsv" ]] && ABSENCES="$FINAL_DIR/absences_explained.tsv"
    SUBTITLE="Final pangenome (after the removal of curated genes)"
else
    echo "ERROR: SET must be preliminary or final" >&2; exit 1
fi
[[ -f "$MATRIX" ]] || { echo "ERROR: missing $MATRIX" >&2; exit 1; }
[[ -d "$GFF_DIR" ]] || { echo "ERROR: missing folder $GFF_DIR" >&2; exit 1; }
echo "set $SET | matrix $MATRIX | GFFs $GFF_DIR | absences ${ABSENCES:-none}"

O="$OUT_DIR/$SET"
mkdir -p "$O/data"
EXTRA=()
[[ -n "$ABSENCES" ]] && EXTRA+=(--absences "$ABSENCES")
[[ -f "$WINNERS" ]] && EXTRA+=(--winners "$WINNERS")
[[ -f "$ENTRY_LIST" ]] && EXTRA+=(--entry-list "$ENTRY_LIST")
python3 "$SD/build_explorer_data.py" --genes-matrix "$MATRIX" --gff-dir "$GFF_DIR" ${EXTRA[@]+"${EXTRA[@]}"} -o "$O/data" \
  --core-missing "$CORE_MISSING" --softcore-missing "$SOFTCORE_MISSING" --cloud-max "$CLOUD_MAX"
python3 "$SD/build_explorer_locus_index.py" --genes-matrix "$MATRIX" --gff-dir "$GFF_DIR" -o "$O/data"

REFERENCE=$(grep -vE '^\s*(#|$)' "$GENOME_ORDER" | head -1)
FATE="$STEP05_DIR/gene_fate_prefer_panel_model.tsv"
if [[ -f "$FATE" && -d "$STEP05_DIR/input/native" ]]; then
    echo "== own annotation layer ($SOURCE_DB)"
    GO=(); [[ -d "$GENOME_DIR" ]] && GO=(--genome-dir "$GENOME_DIR")
    [[ "$SET" == final && -f "$FINAL_DIR/removed_genes.tsv" ]] && GO+=(--removed "$FINAL_DIR/removed_genes.tsv")   # why an original was removed in step 07
    python3 "$SD/build_source_layer.py" --gene-fate "$FATE" --own-gff-dir "$STEP05_DIR/input/native" --gff-dir "$GFF_DIR" \
        --reference "$REFERENCE" ${GO[@]+"${GO[@]}"} -o "$O/data"
else
    echo "own annotation layer: left out (no $FATE or $STEP05_DIR/input/native)"
fi
if [[ -n "$TE_DIR" ]]; then
    [[ -d "$TE_DIR" ]] || { echo "ERROR: no folder $TE_DIR" >&2; exit 1; }
    echo "== TE line ($TE_DIR)"
    python3 "$SD/build_te_layer.py" --te-dir "$TE_DIR" --genomes "$GENOME_ORDER" -o "$O/data"
else
    echo "TE line: left out (TE_DIR not set)"
fi
NG=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))['genomes']))" "$O/data/manifest.json")
CURATION=$([[ "$SET" == preliminary ]] && echo yes || echo no)   # the curation flags: preliminary page only
sed -e "s|__CURATION__|${CURATION}|g" -e "s|__TITLE__|${TITLE}|g" -e "s|__SUBTITLE__|${SUBTITLE}|g" -e "s|__SPECIES__|${SPECIES}|g" \
    -e "s|__NGENOMES__|${NG}|g" -e "s|__SOURCE_DB__|${SOURCE_DB}|g" -e "s|__SOURCE_URL__|${SOURCE_URL}|g" \
    -e "s|__REFERENCE__|${REFERENCE}|g" -e "s|__LINK_GENOMES__|${LINK_GENOMES}|g" "$SD/explorer.html" > "$O/explorer.html"
cp "$SD/serve_explorer.py" "$O/"
if [[ -n "$OG_TABLE" && ! -f "$OG_TABLE" ]]; then
    echo "WARNING: no $OG_TABLE - the page is built without orthogroups" >&2
elif [[ -n "$OG_TABLE" ]]; then
    echo "== orthogroups ($OG_TABLE)"
    rm -rf "$O.og"
    python3 "$SD/add_orthogroups.py" --explorer-dir "$O" --og-table "$OG_TABLE" --reference "$REFERENCE" --out-dir "$O.og"
    rm -rf "$O" && mv "$O.og" "$O"
fi
python3 "$SD/make_standalone_explorer.py" "$O" "$O/explorer_standalone.html"

echo "Done. Copy $O/explorer_standalone.html to your computer and open it with a double click."
