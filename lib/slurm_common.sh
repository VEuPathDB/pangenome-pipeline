#!/usr/bin/env bash
#
# Shared helpers for the pipeline's job scripts. Not run on its own - each job script finds and
# sources it at the top (PIPELINE_DIR, else the checkout the script sits in, else a folder above
# the submit folder):
#
#     _pd="${PIPELINE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
#     [[ -f "$_pd/lib/slurm_common.sh" ]] || { _pd="${SLURM_SUBMIT_DIR:-$PWD}"; while [[ "$_pd" != / && ! -f "$_pd/lib/slurm_common.sh" ]]; do _pd="$(dirname "$_pd")"; done; }
#     source "$_pd/lib/slurm_common.sh" || { echo "ERROR: pipeline not found - export PIPELINE_DIR=/path/to/pangenome-pipeline" >&2; exit 1; }
#
# The pipeline checkout and the species project are separate folders. Jobs are submitted from the
# project folder (genomes/, gff_originals/, genome_order.txt), and every step writes under
# $OUTPUT_DIR/stepNN_<name>/ ($STEP01_OUT ... $STEP08_OUT below).
#
# This file puts the pipeline's Python package (src/) on PYTHONPATH, loads pipeline_config.sh
# (the pipeline's, then the project's if it has one - cluster paths and environment names; copy
# pipeline_config.example.sh to start), sets the output folders, then provides:
#
#   find_script VAR <path/from/pipeline/root>  set VAR to the path of a pipeline file (a shared shell
#                                             script, or a loose script outside src/)
#   activate_env NAME                         activate a micromamba environment
#   require_config VAR [VAR ...]              stop with a clear message if a setting is empty

# Symlinks resolved, so a lib/ linked into a project still leads back to the real checkout.
PIPELINE_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)"
PIPELINE_CONFIG="${PIPELINE_CONFIG:-$PIPELINE_DIR/pipeline_config.sh}"

# The species project: the folder the job was submitted from.
PROJECT_DIR="${PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"

# The pipeline's Python package (src/pangenome_pipeline), for `python3 -m pangenome_pipeline.<step>.<module>`
# with whichever python3 the job has loaded - nothing needs installing.
export PYTHONPATH="$PIPELINE_DIR/src${PYTHONPATH:+:$PYTHONPATH}"

if [[ -f "$PIPELINE_CONFIG" ]]; then
    # shellcheck source=../pipeline_config.example.sh
    source "$PIPELINE_CONFIG"
fi
if [[ -f "$PROJECT_DIR/pipeline_config.sh" && ! "$PROJECT_DIR/pipeline_config.sh" -ef "$PIPELINE_CONFIG" ]]; then
    # shellcheck source=../pipeline_config.example.sh
    source "$PROJECT_DIR/pipeline_config.sh"
fi

# Environment names, used by the wrappers' activate_env calls. pipeline_config.sh can rename them.
: "${ENV_SYNTENY:=synteny}"                 # python3 + pandas/matplotlib/numpy, minimap2, skani
: "${ENV_BIOINFO:=bioinfo}"                 # samtools, seqkit, AGAT
: "${ENV_REPEATMASKER:=repeatmasker4.2}"    # RepeatMasker
: "${ENV_TE_ANNOT:=te_annot}"               # TE-Aid's dependencies (R, BLAST)
: "${ENV_MCHELPER:=mchelper}"               # MCHelper

# Where every step writes: $OUTPUT_DIR/stepNN_<name>/. OUTPUT_DIR defaults to output/ in the project
# and can be set in pipeline_config.sh or per job (sbatch --export=ALL,OUTPUT_DIR=...); a relative
# path is taken from the project folder. Each STEPNN_OUT can also be pointed elsewhere on its own.
: "${OUTPUT_DIR:=output}"
[[ "$OUTPUT_DIR" == /* ]] || OUTPUT_DIR="$PROJECT_DIR/$OUTPUT_DIR"
: "${STEP01_OUT:=$OUTPUT_DIR/step01_genome_selection_qc}"
: "${STEP02_OUT:=$OUTPUT_DIR/step02_te_annotation_clean_up}"
: "${STEP03_OUT:=$OUTPUT_DIR/step03_denovo_genome_annotation}"
: "${STEP04_OUT:=$OUTPUT_DIR/step04_interproscan}"
: "${STEP05_OUT:=$OUTPUT_DIR/step05_cross_genome_annotation_integration}"
: "${STEP06_OUT:=$OUTPUT_DIR/step06_get_pangenes}"
: "${STEP07_OUT:=$OUTPUT_DIR/step07_post_get_pangenes_adjusting}"
: "${STEP08_OUT:=$OUTPUT_DIR/step08_pangenome_visualisation}"
: "${ALT_PANLOCI_OUT:=$OUTPUT_DIR/alternative_get_panloci}"

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
    _pipeline_die "${missing[*]} not set. Set it in $PIPELINE_CONFIG or $PROJECT_DIR/pipeline_config.sh" \
        "(start from $PIPELINE_DIR/pipeline_config.example.sh), or export it before submitting."
}

# find_script VAR REL - REL is a path from the pipeline root, e.g.
# 06_GET_PANGENES/scripts/get_pangenes_common.sh.
find_script() {
    local __var="$1" rel="$2"
    [[ -e "$PIPELINE_DIR/$rel" ]] || _pipeline_die "cannot find $rel in the pipeline ($PIPELINE_DIR)"
    printf -v "$__var" '%s' "$PIPELINE_DIR/$rel"
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
