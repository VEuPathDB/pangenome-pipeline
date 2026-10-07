#!/bin/bash
#SBATCH --job-name=te_kill_list
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=te_kill_list_%j.out
#SBATCH --error=te_kill_list_%j.err
#
# Removes TE-derived gene models from GFFs using a kill list, in one job:
#   1. build_te_kill_lists.py - reads gene_te_overlap.tsv, writes a kill list of
#      gene IDs (TE_overlap_pct >= MIN_PCT) for AGAT.
#   2. AGAT (agat_sp_filter_feature_from_kill_list.pl) - removes those genes from
#      every GFF in GFF_INDIR, writing cleaned copies to GFF_OUTDIR.
#
# Main case: GFF_INDIR="gff_originals" (default below) - clean the original
# per-genome GFFs before running Liftoff.
# Already ran Liftoff? Point GFF_INDIR at "liftoff_polished" instead - the same
# kill list applies, since Liftoff preserves the source gene IDs.
#
# Submit from your project folder:
#     sbatch 02_TE_annotation_Clean_Up/scripts/build_te_kill_lists.slurm.sh
#
# ─── ADAPT: the paths below, the ### ADAPT ### blocks, and the #SBATCH lines ───
GENE_TABLE="te_annotation_integration/gene_te_overlap.tsv"
KILL_OUTDIR="kill_lists"
MIN_PCT=80                        # minimum TE_overlap_pct to kill a protein_coding gene
EXTRA_GENES=""                    # optional: genome_label<TAB>gene_id TSV, no header
GFF_INDIR="gff_originals"         # GFFs to filter - gff_originals (main) or liftoff_polished (alt)
GFF_OUTDIR="gff_cleaned"          # filtered copies land here
# ────────────────────────────────────────────────────────────────────────────

set -euo pipefail

source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1  # find_script, activate_env, pipeline_config.sh
find_script PY 02_TE_annotation_Clean_Up/scripts/build_te_kill_lists.py

[[ -f "$GENE_TABLE" ]] || { echo "ERROR: gene table not found: $GENE_TABLE" >&2; exit 1; }
[[ -d "$GFF_INDIR" ]]  || { echo "ERROR: GFF folder not found: $GFF_INDIR" >&2; exit 1; }
mkdir -p "$GFF_OUTDIR"

### ADAPT ### — put a python3 with pandas on the PATH:
activate_env "$ENV_SYNTENY"  # env names / backend: pipeline_config.sh
### END ADAPT ###

command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }
python3 -c "import pandas" 2>/dev/null || { echo "ERROR: python3 is missing pandas" >&2; exit 1; }

echo "=== 1/2: build_te_kill_lists.py ==="
KILL_ARGS=(--gene-table "$GENE_TABLE" --outdir "$KILL_OUTDIR" --min-pct "$MIN_PCT")
[[ -n "${EXTRA_GENES:-}" ]] && KILL_ARGS+=(--extra-genes "$EXTRA_GENES")
python3 "$PY" "${KILL_ARGS[@]}"
KILL_LIST="$KILL_OUTDIR/te_derived_genes.kill_list.txt"

echo
echo "=== 2/2: AGAT filter (${GFF_INDIR} -> ${GFF_OUTDIR}) ==="
### ADAPT ### — put AGAT (agat_sp_filter_feature_from_kill_list.pl) on the PATH.
#              Re-exports the hook so this activation doesn't depend on the
#              first ADAPT block above still being in scope:
activate_env "$ENV_BIOINFO"  # env names / backend: pipeline_config.sh
### END ADAPT ###

command -v agat_sp_filter_feature_from_kill_list.pl >/dev/null \
    || { echo "ERROR: agat_sp_filter_feature_from_kill_list.pl not on PATH" >&2; exit 1; }

n=0
for gff in "$GFF_INDIR"/*.gff "$GFF_INDIR"/*.gff3; do
    [[ -e "$gff" ]] || continue
    agat_sp_filter_feature_from_kill_list.pl \
        --gff "$gff" \
        --kill_list "$KILL_LIST" \
        --output "$GFF_OUTDIR/$(basename "$gff")"
    n=$((n + 1))
done
(( n > 0 )) || { echo "ERROR: no *.gff/*.gff3 files in $GFF_INDIR" >&2; exit 1; }

echo
echo "Done — kill list in $KILL_OUTDIR/, $n filtered GFF(s) in $GFF_OUTDIR/"
