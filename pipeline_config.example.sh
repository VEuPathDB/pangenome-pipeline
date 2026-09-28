#!/usr/bin/env bash
#
# Pipeline configuration: where the tools, databases and environments live on your cluster.
#
# Copy this file to pipeline_config.sh in the same folder (the project root) and fill it in:
#
#     cp pipeline_config.example.sh pipeline_config.sh
#
# pipeline_config.sh is not committed (see .gitignore). Every job script loads it through
# lib/slurm_common.sh. Anything set here can also be overridden for a single job:
#
#     sbatch --export=ALL,LABEL=Toxo7a 08_Pangenome_Visualisation/scripts/plot_all_job.sh
#
# Only the settings a job actually uses need to be filled in; a job stops with a message
# naming any setting it needs that is still empty.

# ── micromamba ──────────────────────────────────────────────────────────────────────────────
# e.g. on our cluster: "$HOME/hc-storage/micromamba" and "$HOME/hc-storage/bin/micromamba/micromamba"
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
