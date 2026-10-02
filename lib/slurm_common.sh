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
#   activate_env NAME                         put environment NAME's tools on PATH (container,
#                                             micromamba, conda or what's already there - see below)
#   require_config VAR [VAR ...]              stop with a clear message if a setting is empty

# Symlinks resolved, so a lib/ linked into a project still leads back to the real checkout.
export PIPELINE_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")/.." && pwd)"   # exported: jobs a job submits find it too
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

# How activate_env provides them (see activate_env below) and the venv for the pipeline's own Python.
: "${ENV_BACKEND:=auto}"
: "${PIPELINE_VENV:=$PIPELINE_DIR/.venv}"

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
# workflow/step06_get_pangenes/get_pangenes_common.sh.
find_script() {
    local __var="$1" rel="$2"
    [[ -e "$PIPELINE_DIR/$rel" ]] || _pipeline_die "cannot find $rel in the pipeline ($PIPELINE_DIR)"
    printf -v "$__var" '%s' "$PIPELINE_DIR/$rel"
}

# activate_env NAME - put environment NAME's tools on PATH. ENV_BACKEND (pipeline_config.sh) picks
# how; the default, auto, decides per environment:
#
#   container   $CONTAINER_DIR/NAME.sif, or a folder $CONTAINER_DIR/NAME/ of .sif images (e.g. one
#               biocontainers image per tool). Each program in the images' /usr/local/bin and
#               /opt/conda/bin gets a small wrapper on PATH that runs it with apptainer/singularity.
#   micromamba  MAMBA_EXE + MAMBA_ROOT_PREFIX are set.
#   conda       conda is on PATH (or CONDA_EXE is set).
#   path        none of the above: use whatever is already on PATH (module load, ~/.bashrc, ...).
#
# Python: if PIPELINE_VENV (default .venv in the pipeline checkout, made with `uv sync --frozen --no-dev`)
# exists, its python3 is used wherever the environment doesn't bring its own. Works whether or not
# the caller runs with `set -u`.
activate_env() {
    local name="$1" backend="${ENV_BACKEND:-auto}" restore_u=""
    [[ $- == *u* ]] && restore_u=1
    set +u
    if [[ "$backend" == auto ]]; then
        if [[ -n "$(_container_images "$name")" ]]; then
            backend=container
        elif [[ -n "${MAMBA_EXE:-}" ]]; then
            backend=micromamba
        elif [[ -n "${CONDA_EXE:-}" ]] || command -v conda >/dev/null 2>&1; then
            backend=conda
        else
            backend=path
        fi
    fi

    # The venv first, so that anything the environment brings comes before it on PATH.
    if [[ -x "${PIPELINE_VENV:-}/bin/python3" ]]; then
        export PATH="$PIPELINE_VENV/bin:$PATH"
    fi

    case "$backend" in
        container) _activate_container "$name" ;;
        micromamba)
            require_config MAMBA_EXE MAMBA_ROOT_PREFIX
            [[ -x "$MAMBA_EXE" ]] || _pipeline_die "MAMBA_EXE is not an executable: $MAMBA_EXE"
            export MAMBA_EXE MAMBA_ROOT_PREFIX
            eval "$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX")"
            micromamba activate "$name" || _pipeline_die "could not activate micromamba environment '$name'"
            ;;
        conda)
            eval "$("${CONDA_EXE:-conda}" shell.bash hook)"
            conda activate "$name" || _pipeline_die "could not activate conda environment '$name'"
            ;;
        path) ;;
        *) _pipeline_die "ENV_BACKEND must be auto, container, micromamba, conda or path (got '$backend')" ;;
    esac
    echo "activate_env $name: $backend" >&2
    [[ -n "$restore_u" ]] && set -u
    return 0
}

# The .sif images for environment NAME, one per line (none if CONTAINER_DIR is unset).
_container_images() {
    [[ -n "${CONTAINER_DIR:-}" ]] || return 0
    if [[ -f "$CONTAINER_DIR/$1.sif" ]]; then
        echo "$CONTAINER_DIR/$1.sif"
    elif [[ -d "$CONTAINER_DIR/$1" ]]; then
        find "$CONTAINER_DIR/$1" -maxdepth 1 -name '*.sif' | sort
    fi
}

# A folder of wrapper scripts, one per program in NAME's images, put first on PATH. The first image
# (alphabetically) providing a program wins. The project, the pipeline, the output folder and
# CONTAINER_BINDS (comma-separated) are bound into the container.
_activate_container() {
    local name="$1" img tool shim_dir binds="" dir
    local -a images
    mapfile -t images < <(_container_images "$name")
    ((${#images[@]})) || _pipeline_die "no container for '$name': expected $CONTAINER_DIR/$name.sif or $CONTAINER_DIR/$name/*.sif"
    if [[ -z "${CONTAINER_CMD:-}" ]]; then
        command -v apptainer >/dev/null 2>&1 || command -v singularity >/dev/null 2>&1 \
            || { module load apptainer || module load singularity; } >/dev/null 2>&1 || true
        if command -v apptainer >/dev/null 2>&1; then CONTAINER_CMD=apptainer
        elif command -v singularity >/dev/null 2>&1; then CONTAINER_CMD=singularity
        else _pipeline_die "neither apptainer nor singularity found (needed for the '$name' container)"
        fi
    fi

    mkdir -p "$OUTPUT_DIR"
    for dir in "$PROJECT_DIR" "$PIPELINE_DIR" "$OUTPUT_DIR"; do
        binds+="${binds:+,}$(readlink -f "$dir")"
    done
    [[ -n "${CONTAINER_BINDS:-}" ]] && binds+=",$CONTAINER_BINDS"

    shim_dir="$(mktemp -d "${TMPDIR:-/tmp}/pipeline_env_${name}.XXXXXX")"
    for img in "${images[@]}"; do
        while read -r tool; do
            [[ "$tool" =~ ^[A-Za-z0-9._+-]+$ && ! -e "$shim_dir/$tool" ]] || continue
            printf '#!/bin/bash\nexec %q --silent exec --bind %q %q %q "$@"\n' \
                "$CONTAINER_CMD" "$binds" "$img" "$tool" >"$shim_dir/$tool"
            chmod +x "$shim_dir/$tool"
        done < <("$CONTAINER_CMD" exec "$img" sh -c \
            'for d in /usr/local/bin /opt/conda/bin; do [ -d "$d" ] && ls "$d"; done' 2>/dev/null)
    done
    export PATH="$shim_dir:$PATH"
}
