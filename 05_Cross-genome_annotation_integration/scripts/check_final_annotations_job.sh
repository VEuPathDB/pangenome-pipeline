#!/bin/bash
#SBATCH --job-name=check_final_annotations
#SBATCH --partition=low
#SBATCH --output=expanded_annotation/logs/%x_%j.out
#SBATCH --error=expanded_annotation/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#
# Checks that every genome has its <genome>_final_annotation.gff after build_final_annotation_basic.sh.
# For each genome that does not, it says how far the build got, which log to read, and prints the
# sbatch command that rebuilds only those genomes.
#
# It only looks at files and takes a second, so it can be run directly on the login node, from
# your project folder:
#     bash 05_Cross-genome_annotation_integration/scripts/check_final_annotations_job.sh
# or submitted as a job (the log folder must exist first):
#     mkdir -p expanded_annotation/logs && sbatch 05_Cross-genome_annotation_integration/scripts/check_final_annotations_job.sh

set -euo pipefail

# ─── Edit these for your run ───
GENOMES="genome_order.txt"                  # file with one genome name per line (same as the build's ORDER)
OUTPUT_DIR="expanded_annotation"            # build_final_annotation_basic.sh's OUTPUT_DIR
BUILD_SCRIPT="05_Cross-genome_annotation_integration/scripts/build_final_annotation_basic.sh"   # only used to print the resubmit command
# ────────────────────────────────

[[ -f "$GENOMES" ]]    || { echo "ERROR: genome list not found: $GENOMES (run from the project folder)" >&2; exit 1; }
[[ -d "$OUTPUT_DIR" ]] || { echo "ERROR: $OUTPUT_DIR not found - has the build been submitted yet?" >&2; exit 1; }

genomes=()
while IFS= read -r line; do
    line="${line%%#*}"; line="${line//[[:space:]]/}"
    [[ -n "$line" ]] && genomes+=("$line")
done < "$GENOMES"

missing_tasks=()
running_tasks=()
n_ok=0
for i in "${!genomes[@]}"; do
    g="${genomes[$i]}"
    task=$((i + 1))
    if [[ -s "$OUTPUT_DIR/${g}_final_annotation.gff" ]]; then
        echo "OK       $g"
        n_ok=$((n_ok + 1))
        continue
    fi

    # How far did it get? Each layer leaves its own file in the genome's folder.
    if   [[ -s "$OUTPUT_DIR/$g/step_rescued.gff" ]]; then where="the final overlap check"
    elif [[ -s "$OUTPUT_DIR/$g/step_good.gff"    ]]; then where="round 2 (rescue)"
    elif [[ -s "$OUTPUT_DIR/$g/step_main.gff"    ]]; then where="round 1 (good transfers)"
    elif [[ -d "$OUTPUT_DIR/$g"                  ]]; then where="the main step"
    else                                                  where=""
    fi

    # The job's own log (found through the line the build prints for this genome).
    joblog=$(grep -l -E "building ${g}\$|=== Genome: ${g} |=== Genome [0-9]+/[0-9]+: ${g} ===" "$OUTPUT_DIR"/logs/build_final_annotation_*.out 2>/dev/null | tail -1 || true)

    # Is that job still queued or running? Then the genome is not missing, just not finished yet.
    running=0
    if [[ -n "$joblog" ]] && command -v squeue >/dev/null 2>&1; then
        jobid=$(basename "$joblog" .out); jobid="${jobid##*_}"
        if squeue -h -u "${USER:-$(id -un)}" -o "%A" 2>/dev/null | grep -qx "$jobid"; then running=1; fi
    fi

    if [[ $running -eq 1 ]]; then
        echo "RUNNING  $g  (array task $task) - still running, currently in ${where:-the start}"
        running_tasks+=("$task")
        continue
    fi

    if [[ -z "$where" ]]; then stage="never started (no folder)"; else stage="stopped in $where"; fi
    echo "MISSING  $g  (array task $task) - $stage"
    # A failure is reported in the job's .err, e.g. "ERROR: missing ..." or the end of AGAT's log.
    if [[ -n "$joblog" && -s "${joblog%.out}.err" ]]; then
        echo "           job log: ${joblog%.out}.err - last lines:"
        tail -3 "${joblog%.out}.err" | sed 's/^/             | /'
    elif [[ -n "$joblog" ]]; then
        echo "           job log: ${joblog%.out}.err (empty)"
    else
        echo "           no job log found for this genome in $OUTPUT_DIR/logs/ - was it submitted?"
    fi
    latest=$(ls -t "$OUTPUT_DIR/$g"/*.log 2>/dev/null | head -1 || true)
    [[ -n "$latest" ]] && echo "           newest AGAT log: $latest"
    missing_tasks+=("$task")
done

echo
echo "${n_ok} of ${#genomes[@]} genomes have their final annotation" \
     "(${#running_tasks[@]} still running, ${#missing_tasks[@]} missing)."

if [[ ${#missing_tasks[@]} -gt 0 ]]; then
    tasks=$(IFS=,; echo "${missing_tasks[*]}")
    echo
    echo "The lines under each MISSING genome come from its job log. For more, read the end of the AGAT log:"
    echo "    tail -20 <newest AGAT log>"
    echo
    echo "To rebuild only the missing ones:"
    echo "    sbatch --array=${tasks} ${BUILD_SCRIPT}"
    exit 1
fi

if [[ ${#running_tasks[@]} -gt 0 ]]; then
    echo "Nothing is missing, but some genomes are still running - run this again when they finish."
    exit 2
fi

echo "Done - every genome is complete."
