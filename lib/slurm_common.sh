#!/usr/bin/env bash
#
# Shared helpers for the pipeline's job scripts. Not run on its own - each wrapper sources it
# from the project folder it is submitted from:
#
#     source "${SLURM_SUBMIT_DIR:-$PWD}/lib/slurm_common.sh" || exit 1
#
# It loads <project>/pipeline_config.sh if there is one (cluster paths and environment names -
# copy pipeline_config.example.sh to start), then provides:
#
#   find_script VAR <path/from/project/script.py>   set VAR to the path of a pipeline script
#   activate_env NAME                         put environment NAME's tools on PATH (container,
#                                             micromamba, conda or what's already there - see below)
#   require_config VAR [VAR ...]              stop with a clear message if a setting is empty

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PIPELINE_CONFIG="${PIPELINE_CONFIG:-$PIPELINE_DIR/pipeline_config.sh}"

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
: "${ENV_TE_FILTER:=te_filter}"             # BLAST+, cd-hit, seqkit (TE library Step 2/4; only as a container)
: "${ENV_GET_PANGENES:=get_pangenes}"       # get_pangenes.pl (unless GET_PANGENES_ENV is set)

# How activate_env provides them (see activate_env below) and the venv for the pipeline's own Python.
: "${ENV_BACKEND:=auto}"
: "${PIPELINE_VENV:=$PIPELINE_DIR/.venv}"

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

# activate_env NAME - put environment NAME's tools on PATH. ENV_BACKEND (pipeline_config.sh) picks
# how; the default, auto, decides per environment:
#
#   container   $CONTAINER_DIR/NAME.sif, or a folder $CONTAINER_DIR/NAME/ of .sif images (e.g. one
#               biocontainers image per tool). Each program on the images' own PATH (system folders
#               such as /usr/bin excepted) gets a small wrapper that runs it with apptainer/singularity.
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

# Run inside an image: "name path" for each executable on its own PATH, minus the system folders.
# path is the real file when name is a symlink to a file of the same name (bioconda's TE-Aid finds
# its R scripts next to $0, so it must be called at its real location), else just name.
_CONTAINER_LIST_TOOLS='IFS=:; for d in $PATH; do
    case "$d" in /bin|/sbin|/usr/bin|/usr/sbin|/usr/local/sbin|"") continue ;; esac
    for f in "$d"/*; do
        [ -f "$f" ] && [ -x "$f" ] || continue
        n="${f##*/}"; r="$(readlink -f "$f")"
        if [ -L "$f" ] && [ "${r##*/}" = "$n" ]; then echo "$n $r"; else echo "$n $n"; fi
    done
done'

# A folder of wrapper scripts, one per program in NAME's images, put first on PATH. The first image
# (alphabetically) providing a program wins; programs the host system has in /usr/bin or /bin (awk,
# gzip, perl, python3, ...) are left to the host. The submit folder (the project), the pipeline and
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

    for dir in "${SLURM_SUBMIT_DIR:-$PWD}" "$PIPELINE_DIR"; do
        binds+="${binds:+,}$(readlink -f "$dir")"
    done
    [[ -n "${CONTAINER_BINDS:-}" ]] && binds+=",$CONTAINER_BINDS"

    shim_dir="$(mktemp -d "${TMPDIR:-/tmp}/pipeline_env_${name}.XXXXXX")"
    for img in "${images[@]}"; do
        while read -r tool target; do
            [[ "$tool" =~ ^[A-Za-z0-9._+-]+$ && ! -e "$shim_dir/$tool" ]] || continue
            [[ -e "/usr/bin/$tool" || -e "/bin/$tool" ]] && continue   # the host's own awk, perl, python3, ...
            printf '#!/bin/bash\nexec %q --silent exec --bind %q %q %q "$@"\n' \
                "$CONTAINER_CMD" "$binds" "$img" "${target:-$tool}" >"$shim_dir/$tool"
            chmod +x "$shim_dir/$tool"
        done < <("$CONTAINER_CMD" exec "$img" sh -c "$_CONTAINER_LIST_TOOLS" 2>/dev/null)
    done
    export PATH="$shim_dir:$PATH"
}
