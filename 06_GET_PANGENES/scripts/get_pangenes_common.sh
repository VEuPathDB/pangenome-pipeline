#!/bin/bash
#
# Shared body of the step-06 jobs. Not submitted on its own: sourced by
#   link_for_get_pangenes_7a.slurm.sh / link_for_get_pangenes_7b.slurm.sh  -> link_for_get_pangenes
#   run_get_pangenes_7a_job.sh        / run_get_pangenes_7b_job.sh         -> run_get_pangenes
# Each wrapper keeps its own #SBATCH header and "Edit these for your run" settings, then calls
# the matching function below.

# Symlink every genome's GFF and FASTA into OUT_DIR under the same basename, which is how
# get_pangenes.pl pairs them up. Nothing is copied or rewritten.
# Needs: ANNOTATION_DIR ANNOTATION_SUFFIX GENOME_DIR GENOME_SUFFIX ORDER OUT_DIR
link_for_get_pangenes() {
  local GENOMES genome gff fna n=0

  if [[ -f "$ORDER" ]]; then
    mapfile -t GENOMES < <(grep -vE '^\s*(#|$)' "$ORDER")
  else
    IFS=',' read -ra GENOMES <<< "$ORDER"
  fi
  [[ ${#GENOMES[@]} -ge 1 ]] || { echo "ERROR: ORDER has no genomes - check it's been edited" >&2; exit 1; }

  mkdir -p "$OUT_DIR"

  for genome in "${GENOMES[@]}"; do
    gff="$ANNOTATION_DIR/${genome}${ANNOTATION_SUFFIX}"
    fna="$GENOME_DIR/${genome}${GENOME_SUFFIX}"
    [[ -f "$gff" ]] || { echo "ERROR: missing $gff" >&2; exit 1; }
    [[ -f "$fna" ]] || { echo "ERROR: missing $fna" >&2; exit 1; }
    ln -sf "$(cd "$(dirname "$gff")" && pwd)/$(basename "$gff")" "$OUT_DIR/${genome}_for_get_pangenes.gff"
    ln -sf "$(cd "$(dirname "$fna")" && pwd)/$(basename "$fna")" "$OUT_DIR/${genome}_for_get_pangenes.fna"
    n=$((n + 1))
    echo "linked $genome"
  done

  echo "$n genome(s) linked into $OUT_DIR. Outputs are in $OUT_DIR"
}

# Run get_pangenes.pl on INPUT_DIR from the get_pangenes micromamba environment.
# Needs: GET_PANGENES_ENV INPUT_DIR N_CPUS TAXA_THRESHOLD VARIANT; optional TAXON_REGEX
run_get_pangenes() {
  # Remove any previously loaded cluster modules.
  module purge

  # Avoid mixing Perl libraries from other environments.
  unset PERL5LIB
  unset PERL_LOCAL_LIB_ROOT
  unset PERL_MB_OPT
  unset PERL_MM_OPT

  # Add the micromamba environment's own programs (get_pangenes.pl, its bundled Python, bedtools) to PATH.
  export PATH="$GET_PANGENES_ENV:$PATH"

  local GETP="$GET_PANGENES_ENV/get_pangenes.pl"
  local BEDTOOLS="$GET_PANGENES_ENV/bedtools"

  if [[ -n "${TAXON_REGEX:-}" ]]; then
    "$GETP" -d "$INPUT_DIR" -t "$TAXA_THRESHOLD" -s "$TAXON_REGEX" -n "$N_CPUS" -B "$BEDTOOLS"
  else
    "$GETP" -d "$INPUT_DIR" -t "$TAXA_THRESHOLD" -n "$N_CPUS" -B "$BEDTOOLS"
  fi

  echo "get_pangenes ($VARIANT) completed. Outputs (pangene_matrix*.tab, etc.) are in $INPUT_DIR"
}
