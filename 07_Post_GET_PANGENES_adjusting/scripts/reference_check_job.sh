#!/bin/bash
#SBATCH --job-name=reference_check
#SBATCH --partition=low
#SBATCH --output=reference_check/logs/%x_%j.out
#SBATCH --error=reference_check/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00

set -euo pipefail

# Step 07, reference gene check (after part 2, before part 4): one page, reference_check.html, with the genes the
# pipeline added to the reference, their evidence, and a choice per gene. The choices are downloaded from the page as
# curation_actions.tsv.
#   1. region_rnaseq.py          RNA-seq over each new gene, from the bigWig files under BW_ROOT (skipped if BW_ROOT="";
#                                with RERUN_RNASEQ=no the table of an earlier run is used again)
#   2. build_reference_check.py  the page; intron junctions from JUNCTIONS (left out if JUNCTIONS=""); protein
#                                evidence from PROTEIN_EVIDENCE (protein_evidence_job.sh's table; left out if ""); a
#                                "BLAST in NCBI" link per gene, with its protein from GENOME_DIR/<genome>.fna
# Submit from your project folder; the log folder must exist BEFORE submitting:
#     mkdir -p reference_check/logs
#     sbatch 07_Post_GET_PANGENES_adjusting/scripts/reference_check_job.sh

# ===== Edit these =====  (each can also be set on the command line: sbatch --export=ALL,NAME=value ...)
GENOME_ORDER="${GENOME_ORDER-genome_order.txt}"                                   # the reference is its first line
GENOME="${GENOME-}"                                                               # or name the reference here
STEP05_DIR="${STEP05_DIR-basket_annotation}"                                      # step 05 results
GFF_DIR="${GFF_DIR-$STEP05_DIR/gff_prefer_panel_model_final}"                     # step 05 final GFFs
FIXED_DIR="${FIXED_DIR-pangenome_fixed}"                                          # step 07 part 2
JUNCTIONS="${JUNCTIONS-}"                                                         # intron junctions GFF (.gff or .gff.gz); "" = none
BW_ROOT="${BW_ROOT-}"                                                             # folder with the RNA-seq bigWigs; "" = none
PATTERN="${PATTERN-unique_*unlogged*.bw}"                                         # bigWig file names to use
RERUN_RNASEQ="${RERUN_RNASEQ-yes}"                                                # no = reuse the RNA-seq table already in OUT_DIR
OUT_DIR="${OUT_DIR-reference_check}"
OUT_PAGE="${OUT_PAGE-reference_check.html}"                                       # the page, in OUT_DIR (another name keeps the old page)
GENOME_DIR="${GENOME_DIR-genomes}"                                                # <genome>.fna: the protein of each new gene, for its BLAST link; "" = no link
BLAST_EXCLUDE="${BLAST_EXCLUDE-}"                                        # the species left out of the BLAST link (NCBI's name and taxid)
PROTEIN_EVIDENCE="${PROTEIN_EVIDENCE-}"                                           # protein_evidence_job.sh's protein_evidence.tsv; "" = none
EXPLORER="${EXPLORER-../pangenome_explorer/preliminary/explorer_standalone.html}"   # the explorer, from OUT_DIR: each pangene ID links to it
# ======================

set +u
module purge
### ADAPT ### - python3 (with pyBigWig when BW_ROOT is set)
export MAMBA_ROOT_PREFIX="$HOME/hc-storage/micromamba"
export MAMBA_EXE="$HOME/hc-storage/bin/micromamba/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
micromamba activate bioinfo
### END ADAPT ###
set -u

# Locate the scripts. Under sbatch the wrapper is copied to a spool dir, so ${BASH_SOURCE[0]} is useless there -
# fall back to the standard layout relative to the directory you submitted from.
SD=""
for cand in "$(dirname "${BASH_SOURCE[0]}")" \
            "${SLURM_SUBMIT_DIR:-$PWD}/07_Post_GET_PANGENES_adjusting/scripts" \
            "07_Post_GET_PANGENES_adjusting/scripts"; do
    [[ -f "$cand/build_reference_check.py" ]] && { SD="$cand"; break; }
done
[[ -n "$SD" ]] || { echo "ERROR: cannot find build_reference_check.py (run from the project folder)" >&2; exit 1; }

[[ -n "$GENOME" ]] || GENOME=$(grep -vE '^\s*(#|$)' "$GENOME_ORDER" | head -1)
GFF="$GFF_DIR/${GENOME}_basket_annotation.gff"
for f in "$GFF" "$FIXED_DIR/pangene_matrix_genes.tr.tab"; do
    [[ -f "$f" ]] || { echo "ERROR: no $f" >&2; exit 1; }
done
mkdir -p "$OUT_DIR"
echo "reference: $GENOME ($GFF)"

OPT=()
[[ -f "$STEP05_DIR/winners.tsv" ]] && OPT+=(--winners "$STEP05_DIR/winners.tsv")
[[ -f "$FIXED_DIR/pangene_lists/pangene_quality.tsv" ]] && OPT+=(--quality "$FIXED_DIR/pangene_lists/pangene_quality.tsv")
if [[ -n "$JUNCTIONS" ]]; then
    [[ -f "$JUNCTIONS" ]] || { echo "ERROR: no $JUNCTIONS" >&2; exit 1; }
    OPT+=(--junctions "$JUNCTIONS"); echo "intron junctions: $JUNCTIONS"
else
    echo "intron junctions: none (left out of the page)"
fi
if [[ -n "$GENOME_DIR" ]]; then
    FA=$(ls "$GENOME_DIR/$GENOME".f*a* 2>/dev/null | head -1 || true)
    if [[ -n "$FA" ]]; then OPT+=(--genome-fasta "$FA" --blast-exclude "$BLAST_EXCLUDE"); echo "BLAST links: $FA (left out: ${BLAST_EXCLUDE:-none})"
    else echo "WARNING: no $GENOME_DIR/$GENOME.fna - no BLAST links"; fi
fi
if [[ -n "$PROTEIN_EVIDENCE" ]]; then
    [[ -f "$PROTEIN_EVIDENCE" ]] || { echo "ERROR: no $PROTEIN_EVIDENCE" >&2; exit 1; }
    OPT+=(--protein-evidence "$PROTEIN_EVIDENCE"); echo "protein evidence: $PROTEIN_EVIDENCE"
fi
RT="$OUT_DIR/${GENOME}_region_rnaseq.tsv"
if [[ -n "$BW_ROOT" && "$RERUN_RNASEQ" == no && -f "$RT" ]]; then
    echo "RNA-seq: reusing $RT"
    OPT+=(--region-rnaseq "$RT")
elif [[ -n "$BW_ROOT" ]]; then
    [[ -d "$BW_ROOT" ]] || { echo "ERROR: no folder $BW_ROOT" >&2; exit 1; }
    python3 -c "import pyBigWig" || { echo "ERROR: pyBigWig missing in this environment (see ADAPT)" >&2; exit 1; }
    echo "== 1. RNA-seq over the new genes ($BW_ROOT, $PATTERN)"
    python3 "$SD/region_rnaseq.py" --gff "$GFF" --bw-root "$BW_ROOT" --pattern "$PATTERN" --out "$RT"
    OPT+=(--region-rnaseq "$RT")
else
    echo "RNA-seq: none (left out of the page)"
fi

echo "== 2. the page"
python3 "$SD/build_reference_check.py" --genome "$GENOME" --gff "$GFF" --matrix "$FIXED_DIR/pangene_matrix_genes.tr.tab" \
    ${OPT[@]+"${OPT[@]}"} --explorer "$EXPLORER" --template "$SD/reference_check_template.html" --out "$OUT_DIR/$OUT_PAGE"
echo "Done: $OUT_DIR/$OUT_PAGE"
