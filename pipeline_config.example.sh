#!/usr/bin/env bash
#
# Pipeline configuration: where the tools, databases and environments live on your cluster.
#
# Copy this file to pipeline_config.sh next to it (the pipeline checkout) and fill it in:
#
#     cp pipeline_config.example.sh pipeline_config.sh
#
# A species project can also have its own pipeline_config.sh (in the folder you submit jobs from);
# it is loaded after this one, so it only needs the settings that differ for that project.
# pipeline_config.sh is not committed (see .gitignore). Every job script loads it through
# lib/slurm_common.sh. Anything set here can also be overridden for a single job:
#
#     sbatch --export=ALL,LABEL=Toxo7a "$PIPELINE_DIR"/workflow/step08_pangenome_visualisation/plot_all_job.sh
#
# Only the settings a job actually uses need to be filled in; a job stops with a message
# naming any setting it needs that is still empty.
#
# SLURM partition/account: the #SBATCH lines are read before this file, so they can't come from
# here. Set them for every job in your shell instead (they override the #SBATCH lines):
#     export SBATCH_PARTITION=<partition>  SBATCH_ACCOUNT=<account>

# ── outputs ─────────────────────────────────────────────────────────────────────────────────
# Every step writes under $OUTPUT_DIR/stepNN_<name>/. Default: output/ in the project folder;
# a relative path is taken from the project folder.
OUTPUT_DIR="${OUTPUT_DIR:-output}"

# ── tool environments ───────────────────────────────────────────────────────────────────
# Jobs call `activate_env <name>` (names below). ENV_BACKEND=auto picks, per environment, the first
# of these that is set up - so micromamba is optional:
#   1. container   $CONTAINER_DIR/<name>.sif, or $CONTAINER_DIR/<name>/*.sif (one image per tool is
#                  fine, e.g. biocontainers minimap2 + skani for "synteny"); run with apptainer/singularity
#   2. micromamba  MAMBA_EXE and MAMBA_ROOT_PREFIX below
#   3. conda       conda on PATH
#   4. path        whatever is already on PATH (module load in ~/.bashrc, ...)
# Force one with ENV_BACKEND=container|micromamba|conda|path.
ENV_BACKEND="${ENV_BACKEND:-auto}"
CONTAINER_DIR="${CONTAINER_DIR:-}"
# Extra folders the containers need to see (comma-separated), e.g. where genome symlinks point.
# The project, the pipeline and OUTPUT_DIR are always bound.
CONTAINER_BINDS="${CONTAINER_BINDS:-}"

# The pipeline's own Python code (pandas, matplotlib, ...) runs from this venv when the environment
# doesn't bring a python3. Make it once, pinned to uv.lock, from the pipeline checkout:
#     uv sync --frozen --no-dev
PIPELINE_VENV="${PIPELINE_VENV:-}"   # empty = .venv in the pipeline checkout

# micromamba (optional), e.g. on our cluster: "$HOME/hc-storage/micromamba" and "$HOME/hc-storage/bin/micromamba/micromamba"
MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-}"
MAMBA_EXE="${MAMBA_EXE:-}"

# Environment names (these are the defaults - change them if yours are called something else).
ENV_SYNTENY="${ENV_SYNTENY:-synteny}"                 # python3 + pandas/matplotlib/numpy, minimap2, skani
ENV_BIOINFO="${ENV_BIOINFO:-bioinfo}"                 # samtools, seqkit, AGAT
ENV_REPEATMASKER="${ENV_REPEATMASKER:-repeatmasker4.2}"
ENV_TE_ANNOT="${ENV_TE_ANNOT:-te_annot}"
ENV_MCHELPER="${ENV_MCHELPER:-mchelper}"

# ── step 02: TE annotation ──────────────────────────────────────────────────────────────────
# e.g. "$HOME/hc-storage/bin/TE-Aid/TE-Aid"
TE_AID="${TE_AID:-}"
# e.g. "$HOME/hc-storage/bin/MCHelper/MCHelper.py"
MCHELPER_PY="${MCHELPER_PY:-}"
# e.g. "$HOME/hc-storage/db_TEs/RepBase31.06.fasta/RepBase31.06.fasta"
REPBASE_NT_FASTA="${REPBASE_NT_FASTA:-}"
# e.g. "$HOME/hc-storage/db_TEs/RepeatPeps.lib"
TE_PEP_FASTA="${TE_PEP_FASTA:-}"
# BUSCO lineage HMMs for MCHelper, e.g. ".../busco_downloads/eurotiales_odb12.hmm" - species-specific
BUSCO_LINEAGE="${BUSCO_LINEAGE:-}"

# ── step 03: funannotate ────────────────────────────────────────────────────────────────────
# Shared funannotate install (containers/, databases/, overlay/),
# e.g. "/mnt/hc-storage/groups/veupathdb_annotation/eb_tools/funannotate"
FUNANNOTATE_SHARED="${FUNANNOTATE_SHARED:-}"

# ── step 06: get_pangenes ───────────────────────────────────────────────────────────────────
# bin/ of the get_pangenes micromamba environment, e.g. "$HOME/hc-storage/micromamba/envs/get_pangenes/bin"
GET_PANGENES_ENV="${GET_PANGENES_ENV:-}"

# ── step 08: plots ──────────────────────────────────────────────────────────────────────────
# InterPro entry.list (IPR id -> name), shared by all projects; downloaded there if missing,
# e.g. "$HOME/hc-storage/interpro_entry.list"
INTERPRO_ENTRY_LIST="${INTERPRO_ENTRY_LIST:-}"
