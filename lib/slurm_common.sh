#!/usr/bin/env bash
#
# Shared helpers for the pipeline's job scripts. Not run on its own - each wrapper sources it
# from the project folder it is submitted from:
#
#     source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1
#
# It puts the pipeline's Python package (src/) on PYTHONPATH, loads <project>/pipeline_config.sh if
# there is one (cluster paths and environment names - copy pipeline_config.example.sh to start),
# then provides:
#
#   find_script VAR <path/from/project/script.py>   set VAR to the path of a pipeline script
#   activate_env NAME                         activate a micromamba environment
#   require_config VAR [VAR ...]              stop with a clear message if a setting is empty

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIPELINE_CONFIG="${PIPELINE_CONFIG:-$PIPELINE_DIR/pipeline_config.sh}"

# The pipeline's Python package (src/pangenome_pipeline), for `python3 -m pangenome_pipeline.<step>.<module>`
# with whichever python3 the job has loaded - nothing needs installing.
export PYTHONPATH="$PIPELINE_DIR/src${PYTHONPATH:+:$PYTHONPATH}"

if [[ -f "$PIPELINE_CONFIG" ]]; then
    # shellcheck source=../pipeline_config.example.sh
    source "$PIPELINE_CONFIG"
fi

# Environment names, used by the wrappers' activate_env calls. pipeline_config.sh can rename them.
: "${ENV_SYNTENY:=synteny}"                 # python3 + pandas/matplotlib/numpy, minimap2, skani
: "${ENV_BIOINFO:=bioinfo}"                 # samtools, seqkit, AGAT
: "${ENV_REPEATMASKER:=repeatmasker4.2}"    # RepeatMasker
: "${ENV_TE_ANNOT:=te_annot}"               # TE-Aid's dependencies (R, BLAST)
: "${ENV_MCHELPER:=mchelper}"               # MCHelper

_pipeline_die() {
    echo "ERROR: $*" >&2
    exit 1
}

# require_config VAR [VAR ...] - every named setting must be non-empty.
require_config() {
    local var missing=()
    for var in "$@"; do
        [[ -n "${!var:-}" ]] || missing+=("$var")
    done
    ((${#missing[@]} == 0)) && return 0
    _pipeline_die "${missing[*]} not set. Set it in $PIPELINE_CONFIG (start from" \
        "$PIPELINE_DIR/pipeline_config.example.sh), or export it before submitting."
}

# find_script VAR REL - REL is a path from the project folder, e.g.
# 05_Cross-genome_annotation_integration/scripts/filter_final_annotation.py. Tried next to the
# running wrapper first, then under the submit folder (under sbatch the wrapper itself runs from
# a spool copy, so its own location is no help), then relative to the current folder.
find_script() {
    local __var="$1" rel="$2" cand
    local wrapper="${BASH_SOURCE[${#BASH_SOURCE[@]} - 1]}"
    for cand in \
        "$(dirname "$wrapper")/$(basename "$rel")" \
        "${SLURM_SUBMIT_DIR:-$PWD}/$rel" \
        "$rel"; do
        if [[ -e "$cand" ]]; then
            printf -v "$__var" '%s' "$cand"
            return 0
        fi
    done
    _pipeline_die "cannot find $(basename "$rel") (submit from the project folder)"
}

# activate_env NAME - micromamba environment NAME, via MAMBA_EXE / MAMBA_ROOT_PREFIX from the
# config. Works whether or not the caller runs with `set -u`.
activate_env() {
    local name="$1" restore_u=""
    require_config MAMBA_EXE MAMBA_ROOT_PREFIX
    [[ -x "$MAMBA_EXE" ]] || _pipeline_die "MAMBA_EXE is not an executable: $MAMBA_EXE"
    export MAMBA_EXE MAMBA_ROOT_PREFIX
    [[ $- == *u* ]] && restore_u=1
    set +u
    eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
    micromamba activate "$name" || _pipeline_die "could not activate micromamba environment '$name'"
    [[ -n "$restore_u" ]] && set -u
    return 0
}
