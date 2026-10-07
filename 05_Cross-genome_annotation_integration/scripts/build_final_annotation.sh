#!/usr/bin/env bash
#
# build_final_annotation.sh
#
# Builds one complete annotation per genome, each anchored on that genome's
# own coordinates, using every other genome in the set as evidence - staged
# into two rounds so a low-confidence transfer never blocks a better one
# from another genome. Needs an all-against-all Liftoff run done already
# (every genome's own annotation lifted onto every other genome's FASTA).
#
# For each genome in turn as TARGET (first genome in -r ORDER = REFERENCE):
#   - if TARGET is the reference: its own annotation is used as-is, tagged
#     "main" (nothing was ever lifted onto it from itself)
#   - otherwise: the reference's annotation lifted onto TARGET is filtered
#     to good transfers (valid ORF, coverage & identity >= threshold,
#     primary copy) -> "main"; complemented with TARGET's own original
#     annotation (unfiltered, tagged with its own source)
#   then, for every OTHER genome (not reference, not target):
#     round 1 (good): that genome's own annotation lifted onto TARGET,
#       filtered to good transfers, tagged "good_<genome>", in priority order
#     round 2 (rescue), only once round 1 is done for every genome: the
#       reference's full unfiltered lift ("rescued_<reference>") first, then
#       every other genome's, unfiltered ("rescued_<genome>"), same order
#
# Splitting into two rounds means a locus only falls back to a low-confidence
# transfer if NO genome could supply a good one - a low-confidence gene from
# one genome can no longer block a better version of the same gene coming
# from another.
#
# Requires: python3, AGAT (agat_sp_filter_feature_from_keep_list.pl,
# agat_sp_complement_annotations.pl) on PATH.
#
# Usage: ./build_final_annotation.sh -i /path/to/genome_files -p /path/to/pairwise_liftoff \
#          -o /path/to/output -r genome1,genome2,genome3
#
set -euo pipefail

INPUT_DIR="."
PAIRWISE_DIR="./pairwise_liftoff"
OUTPUT_DIR="./expanded_annotation"
ORDER=""
MIN_COVERAGE="0.9"
MIN_IDENTITY="0.9"
ORIGINAL_SUFFIX=""
PAIR_SUFFIX=""
TARGET_ONLY=""

usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT_DIR -p PAIRWISE_DIR -o OUTPUT_DIR -r ORDER -g ORIGINAL_SUFFIX -x PAIR_SUFFIX [options]

  -i INPUT_DIR       Directory with "<genome><ORIGINAL_SUFFIX>" - each genome's own
                      annotation, including the reference's (default: current directory)
  -p PAIRWISE_DIR     Directory with "<genome>_onto_<other><PAIR_SUFFIX>" files from
                      an all-against-all Liftoff run (default: ./pairwise_liftoff)
  -o OUTPUT_DIR       Where to write outputs (default: ./expanded_annotation)
  -r ORDER            Genome order, first = reference: either a comma-separated list
                      ("genome1,genome2,genome3") or a path to a file with one
                      genome name per line (required)
  -c MIN_COVERAGE     Minimum Liftoff coverage for a "good" transfer (default: 0.9)
  -d MIN_IDENTITY     Minimum Liftoff sequence identity for a "good" transfer (default: 0.9)
  -g ORIGINAL_SUFFIX  Filename suffix for each genome's own original GFF (required,
                      e.g. "_genomic.gff" - no default, this varies by data source)
  -x PAIR_SUFFIX      Filename suffix for all-against-all Liftoff outputs (required,
                      e.g. "_liftoff.gff3_polished" - no default, depends on how you ran Liftoff)
  -t TARGET_ONLY      Build only this one genome instead of every genome in -r ORDER
                      (must still be listed in ORDER, so the reference is known) - lets
                      you run every genome as a separate parallel job, e.g. a SLURM array
                      (build_final_annotation_basic.sh does this with --array). Omit to build every
                      genome in ORDER, one after another, like a normal single-machine run.
  -h                  Show this help
EOF
  exit 1
}

while getopts "i:p:o:r:c:d:g:x:t:h" opt; do
  case "$opt" in
    i) INPUT_DIR="$OPTARG" ;;
    p) PAIRWISE_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    r) ORDER="$OPTARG" ;;
    c) MIN_COVERAGE="$OPTARG" ;;
    d) MIN_IDENTITY="$OPTARG" ;;
    g) ORIGINAL_SUFFIX="$OPTARG" ;;
    x) PAIR_SUFFIX="$OPTARG" ;;
    t) TARGET_ONLY="$OPTARG" ;;
    h|*) usage ;;
  esac
done

[[ -z "$ORDER" ]] && { echo "ERROR: -r ORDER is required" >&2; usage; }
[[ -z "$ORIGINAL_SUFFIX" ]] && { echo "ERROR: -g ORIGINAL_SUFFIX is required" >&2; usage; }
[[ -z "$PAIR_SUFFIX" ]] && { echo "ERROR: -x PAIR_SUFFIX is required" >&2; usage; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for cmd in python3 agat_sp_filter_feature_from_keep_list.pl agat_sp_complement_annotations.pl; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: '$cmd' not found on PATH. Is AGAT installed/activated (e.g. 'conda activate agat')?" >&2
    exit 1
  fi
done

if [[ -f "$ORDER" ]]; then
  mapfile -t GENOMES < <(grep -vE '^\s*(#|$)' "$ORDER")
else
  IFS=',' read -ra GENOMES <<< "$ORDER"
fi
[[ ${#GENOMES[@]} -ge 2 ]] || { echo "ERROR: need at least 2 genomes in ORDER" >&2; exit 1; }

REFERENCE="${GENOMES[0]}"
mkdir -p "$OUTPUT_DIR"

# Absolute paths: every AGAT call below runs inside the genome's own folder (see agat_run).
for d in INPUT_DIR PAIRWISE_DIR OUTPUT_DIR; do
  [[ -d "${!d}" ]] || { echo "ERROR: $d not found: ${!d}" >&2; exit 1; }
  printf -v "$d" '%s' "$(cd "${!d}" && pwd)"
done

# AGAT writes agat_log_* / agat_tmp_* folders into its working directory, and names the temporary
# one after the reference file plus EVERY --add file joined together - with many genomes and long
# names that passes the 255-character file-name limit ("File name too long"). So:
#  - every AGAT call runs inside the genome's own folder (nothing lands in the project folder, and
#    genomes built in parallel never share a log/temp folder);
#  - files given to --add get short names (short_link) - the real files keep their descriptive names.
agat_run() {  # $1=directory to run in, $2=log file for the AGAT output, then the command
  local dir="$1" log="$2"; shift 2
  local attempt max_attempts=3 out="" prev="" arg
  # the file the command is going to write (its --out / --output argument)
  for arg in "$@"; do
    if [[ "$prev" == "--out" || "$prev" == "--output" ]]; then out="$arg"; fi
    prev="$arg"
  done
  # AGAT occasionally crashes at the end of reading a file ("progress bar already finished at
  # .../AGAT/OmniscientI.pm") - a random glitch in its progress bar, not a problem with the data.
  # A second try normally goes through, so a failed call is repeated before giving up. A call that
  # "succeeds" but leaves no annotation behind is treated as a failure too - accepting an empty
  # result would silently produce an empty annotation.
  # For some inputs the crash repeats on every try, so AGAT's progress bar is also switched off:
  # AGAT reads an agat_config.yaml from the folder it runs in, created here if it is missing.
  if [[ ! -f "$dir/agat_config.yaml" ]]; then
    ( cd "$dir" && agat config --expose >/dev/null 2>&1 ) || true
  fi
  if [[ -f "$dir/agat_config.yaml" ]]; then
    sed -i 's/^\( *progress_bar:\).*/\1 false/' "$dir/agat_config.yaml"
  fi
  for attempt in $(seq 1 "$max_attempts"); do
    # start every attempt clean: no half-written output, temp folder or log folder from the last one
    rm -rf "$dir"/agat_tmp_* 2>/dev/null || true
    if [[ -n "$out" ]]; then
      rm -f "$out" "${out%.gff}_report.txt"
      rm -rf "$dir/agat_log_$(basename "${out%.gff}")" 2>/dev/null || true
    fi
    local ok=0
    if ( cd "$dir" && "$@" ) > "$log" 2>&1; then
      if [[ -z "$out" ]] || { [[ -s "$out" ]] && grep -qv '^#' "$out"; }; then
        return 0
      fi
      ok=1
    fi
    if [[ "$attempt" -lt "$max_attempts" ]]; then
      if [[ $ok -eq 1 ]]; then
        echo "WARNING: $1 finished but wrote no annotation to $(basename "$out") (attempt $attempt of $max_attempts), trying again." >&2
      else
        echo "WARNING: $1 failed (attempt $attempt of $max_attempts), trying again. Last line of its log: $(tail -1 "$log")" >&2
      fi
    fi
  done
  # AGAT's own messages go to its log file, so a failure would otherwise be silent in the job log
  echo "ERROR: $1 failed $max_attempts times - last lines of $log:" >&2
  tail -25 "$log" >&2
  return 1
}
short_link() {  # $1=file $2=short name (same folder): a hard link, or a copy if links aren't allowed
  rm -f "$2"
  ln "$1" "$2" 2>/dev/null || cp "$1" "$2"
}

if [[ -n "$TARGET_ONLY" ]]; then
  found=0
  for g in "${GENOMES[@]}"; do [[ "$g" == "$TARGET_ONLY" ]] && found=1; done
  [[ "$found" -eq 1 ]] || { echo "ERROR: -t $TARGET_ONLY is not one of the genomes in -r ORDER" >&2; exit 1; }
  TARGETS=("$TARGET_ONLY")
else
  TARGETS=("${GENOMES[@]}")
fi

relabel_source() {  # $1=file $2=new_source
  awk -F'\t' -v src="$2" 'BEGIN{OFS="\t"} /^#/{print; next} {$2=src; print}' "$1" > "$1.tmp"
  mv "$1.tmp" "$1"
}

relabel_gene_type() {  # $1=file $2=target_type -> stdout redirected by caller
  awk -F'\t' -v target_type="$2" 'BEGIN{OFS="\t"} /^#/{print; next}
    { if (($3=="gene" || $3=="protein_coding_gene") && $9 ~ /gene_biotype=protein_coding/) $3=target_type; print }' "$1"
}

# Round 1: filter every other genome's pairwise Liftoff file to good transfers, tag
# "good_<other>", then complement them ALL in with a single AGAT call (one --add per
# genome, in priority order - AGAT cascades multiple --add files internally, so this
# gives the exact same result as one-genome-at-a-time, just without paying the cost of
# starting AGAT and re-parsing the growing accumulated file once per genome).
run_good_round() {  # $1=target $2=target_dir -> updates $accumulated
  local target="$1" target_dir="$2"
  local add_args=()

  for other in "${others[@]}"; do
    local pair_gff="$PAIRWISE_DIR/${other}_onto_${target}${PAIR_SUFFIX}"
    if [[ ! -f "$pair_gff" ]]; then
      echo "  WARNING: missing $pair_gff, skipping good_$other" >&2
      continue
    fi
    local keep_ids="$target_dir/${other}_good_ids.txt"
    local filtered="$target_dir/${other}_good.gff"
    python3 "$SCRIPT_DIR/get_good_transfer_ids.py" \
      --gff "$pair_gff" --out "$keep_ids" \
      --min-coverage "$MIN_COVERAGE" --min-identity "$MIN_IDENTITY"
    rm -f "$filtered" "${filtered%.gff}_report.txt"
    agat_run "$target_dir" "$target_dir/${other}_good_filter.log" agat_sp_filter_feature_from_keep_list.pl --gff "$pair_gff" --keep_list "$keep_ids" --output "$filtered"
    relabel_gene_type "$filtered" "$main_gene_type" > "$filtered.tmp"
    mv "$filtered.tmp" "$filtered"
    relabel_source "$filtered" "good_${other}"
    local short="$target_dir/a$(printf '%02d' $(( ${#add_args[@]} / 2 + 1 ))).gff"
    short_link "$filtered" "$short"
    add_args+=(--add "$short")
  done

  [[ ${#add_args[@]} -eq 0 ]] && return

  local next="$target_dir/step_good.gff"
  agat_run "$target_dir" "$target_dir/good_complement.log" agat_sp_complement_annotations.pl --ref "$accumulated" "${add_args[@]}" --out "$next"
  echo "  round 1: $(grep '^We added' "$target_dir/good_complement.log" | tr '\n' ' ')"
  accumulated="$next"
}

# Round 2: relabel every other genome's FULL, unfiltered pairwise Liftoff file (the
# reference's own first), then complement them ALL in with a single AGAT call, same
# priority order as round 1.
run_rescue_round() {  # $1=target $2=target_dir -> updates $accumulated
  local target="$1" target_dir="$2"
  local add_args=()
  local rescue_order=()

  [[ "$target" != "$REFERENCE" ]] && rescue_order+=("$REFERENCE")
  rescue_order+=("${others[@]}")

  for other in "${rescue_order[@]}"; do
    local pair_gff="$PAIRWISE_DIR/${other}_onto_${target}${PAIR_SUFFIX}"
    if [[ ! -f "$pair_gff" ]]; then
      echo "  WARNING: missing $pair_gff, skipping rescued_$other" >&2
      continue
    fi
    local relabeled="$target_dir/${other}_rescue.gff"
    relabel_gene_type "$pair_gff" "$main_gene_type" > "$relabeled"
    relabel_source "$relabeled" "rescued_${other}"
    local short="$target_dir/r$(printf '%02d' $(( ${#add_args[@]} / 2 + 1 ))).gff"
    short_link "$relabeled" "$short"
    add_args+=(--add "$short")
  done

  [[ ${#add_args[@]} -eq 0 ]] && return

  local next="$target_dir/step_rescued.gff"
  agat_run "$target_dir" "$target_dir/rescue_complement.log" agat_sp_complement_annotations.pl --ref "$accumulated" "${add_args[@]}" --out "$next"
  echo "  round 2: $(grep '^We added' "$target_dir/rescue_complement.log" | tr '\n' ' ')"
  accumulated="$next"
}

total=${#TARGETS[@]}
n=0
for TARGET in "${TARGETS[@]}"; do
  n=$((n + 1))
  if [[ -n "$TARGET_ONLY" ]]; then
    echo "=== Genome: $TARGET (single-target mode) ==="
  else
    echo "=== Genome $n/$total: $TARGET ==="
  fi
  target_dir="$OUTPUT_DIR/$TARGET"
  mkdir -p "$target_dir"

  original_gff="$INPUT_DIR/${TARGET}${ORIGINAL_SUFFIX}"
  [[ -f "$original_gff" ]] || { echo "ERROR: missing $original_gff" >&2; exit 1; }

  accumulated="$target_dir/step_main.gff"

  if [[ "$TARGET" == "$REFERENCE" ]]; then
    # The reference's own annotation is already the trusted, complete set - use as-is.
    cp "$original_gff" "$accumulated"
    relabel_source "$accumulated" "main"
    main_gene_type=$(cut -f3 "$accumulated" | grep -E '^(gene|protein_coding_gene)$' | sort -u | head -1)
    [[ -n "$main_gene_type" ]] || { echo "ERROR: could not detect gene feature type for $TARGET" >&2; exit 1; }
  else
    ref_pair="$PAIRWISE_DIR/${REFERENCE}_onto_${TARGET}${PAIR_SUFFIX}"
    [[ -f "$ref_pair" ]] || { echo "ERROR: missing $ref_pair" >&2; exit 1; }

    keep_ids="$target_dir/main_good_ids.txt"
    main_gff="$target_dir/main.gff"
    python3 "$SCRIPT_DIR/get_good_transfer_ids.py" \
      --gff "$ref_pair" --out "$keep_ids" \
      --min-coverage "$MIN_COVERAGE" --min-identity "$MIN_IDENTITY"
    rm -f "$main_gff" "${main_gff%.gff}_report.txt"
    agat_run "$target_dir" "$target_dir/main_filter.log" agat_sp_filter_feature_from_keep_list.pl --gff "$ref_pair" --keep_list "$keep_ids" --output "$main_gff"
    relabel_source "$main_gff" "main"

    main_gene_type=$(cut -f3 "$main_gff" | grep -E '^(gene|protein_coding_gene)$' | sort -u | head -1)
    [[ -n "$main_gene_type" ]] || { echo "ERROR: could not detect gene feature type for $TARGET" >&2; exit 1; }

    original_relabeled="$target_dir/original_relabeled.gff"
    relabel_gene_type "$original_gff" "$main_gene_type" > "$original_relabeled"

    agat_run "$target_dir" "$target_dir/main_complement.log" agat_sp_complement_annotations.pl --ref "$main_gff" --add "$original_relabeled" --out "$accumulated"
    echo "  + original annotation: $(grep '^We added' "$target_dir/main_complement.log" | tr '\n' ' ')"
  fi
  echo "  main done: $accumulated"

  others=()
  for g in "${GENOMES[@]}"; do
    [[ "$g" == "$TARGET" || "$g" == "$REFERENCE" ]] && continue
    others+=("$g")
  done

  echo "  round 1 (good transfers only):"
  run_good_round "$TARGET" "$target_dir"

  echo "  round 2 (rescue, unfiltered):"
  run_rescue_round "$TARGET" "$target_dir"

  final="$OUTPUT_DIR/${TARGET}_final_annotation.gff"
  cp "$accumulated" "$final"
  [[ -s "$final" ]] && grep -qv '^#' "$final" || { echo "ERROR: $final is empty - the last AGAT step left no annotation" >&2; exit 1; }

  echo "  done: $final"
done

echo ""
if [[ -n "$TARGET_ONLY" ]]; then
  echo "Done: $OUTPUT_DIR/${TARGET_ONLY}_final_annotation.gff"
else
  echo "All genomes processed. Final annotations: $OUTPUT_DIR/<genome>_final_annotation.gff"
fi
