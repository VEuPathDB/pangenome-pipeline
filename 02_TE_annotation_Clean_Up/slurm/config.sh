#!/usr/bin/env bash
# Configuration file for the TE candidate discovery pipeline.
#

VERSION="1.1"

###############################################################################
# Project directory
###############################################################################

# The pipeline uses one working directory for both input files and output files.
# Run the pipeline from the directory that contains:
#   - genome FASTA files
#   - curated TE library
#   - db/ folder with RepBase and TE peptide files
#   - pipeline scripts
#
# Output folders will also be created in this same directory.

# This file lives at  <project>/02_TE_annotation_Clean_Up/slurm/config.sh .
# Derive the directories from its own location, so every pipeline script gets the
# same paths no matter which directory it is run from.
_CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"   # .../02_TE_annotation_Clean_Up/slurm
TE_STEP_DIR="$(cd "${_CONFIG_DIR}/.." && pwd -P)"                # .../02_TE_annotation_Clean_Up
PROJECT_DIR="$(cd "${TE_STEP_DIR}/.." && pwd -P)"               # the project root (holds genomes/)
unset _CONFIG_DIR

# Cluster paths (MCHelper, RepBase, TE peptides, BUSCO lineage, micromamba) come from the
# project's pipeline_config.sh - see pipeline_config.example.sh. This also provides activate_env.
source "${PROJECT_DIR}/lib/slurm_common.sh"

###############################################################################
# Genome inputs
###############################################################################
#
# Default: genome FASTA files are expected in the same directory where you run
# the pipeline.
#
# If your genome FASTA files are stored in another directory, replace
# "${PROJECT_DIR}" with the full path to that directory.

GENOME_DIR="${PROJECT_DIR}/genomes"

# Edit GENOME_GLOB and GENOME_SUFFIX only if your genome files use a different
# naming pattern, for example "genome.fasta"
#
# For these names:
#   GENOME_GLOB finds all genome FASTA files.
#   GENOME_SUFFIX is removed from the file name to create the genome ID.
#   
GENOME_GLOB="*.fna" 
GENOME_SUFFIX=".fna"

###############################################################################
# TE libraries and MCHelper inputs
###############################################################################

# Curated TE library used in Step 2 to remove TE families that are already known.
# This can be obtained from Repbase depending on the species. In this case, LTRs should be
# joined to the internal region.
# Edit the file name for your curated library.
CURATED_TE_LIBRARY="${TE_STEP_DIR}/fake_lib.fasta"

# BUSCO lineage used by MCHelper.
# Check the best BUSCO database to your species
#Example on how to obtain the BUSCO dataset for Drosophila
#wget https://busco-data.ezlab.org/v4/data/lineages/diptera_odb10.2020-08-05.tar.gz
#tar xvf diptera_odb10.2020-08-05.tar.gz
#cat diptera_odb10/hmms/*.hmm > diptera_odb10.hmm
#cd -

# Set BUSCO_LINEAGE in pipeline_config.sh.
BUSCO_LINEAGE="${BUSCO_LINEAGE:-}"

# MCHelper script - set MCHELPER_PY in pipeline_config.sh.
MCHELPER_PY="${MCHELPER_PY:-}"

###############################################################################
# Reference files for final BLAST searches
###############################################################################

# RepBase nucleotide library used by BLASTN in the final filtering step.
# Set REPBASE_NT_FASTA in pipeline_config.sh.
REPBASE_NT_FASTA="${REPBASE_NT_FASTA:-}"

# TE peptide library, for example RepeatPeps.lib, used by BLASTX in the final
# filtering step. Set TE_PEP_FASTA in pipeline_config.sh.
TE_PEP_FASTA="${TE_PEP_FASTA:-}"

###############################################################################
# Output directories
###############################################################################

CONFIG_DIR="${TE_STEP_DIR}/config"
LOG_DIR="${TE_STEP_DIR}/logs"

GENOME_TABLE="${CONFIG_DIR}/genomes.tsv"
MCHELPER_TABLE="${CONFIG_DIR}/mchelper_inputs.tsv"

REPEATMODELER_OUTDIR="${TE_STEP_DIR}/repeatmodeler_more_genomes"
CANDIDATE_OUTDIR="${TE_STEP_DIR}/TE_candidate_screen_original_headers"
MCHELPER_OUTDIR="${TE_STEP_DIR}/MCHelper_newTEs"
FINAL_OUTDIR="${TE_STEP_DIR}/final_TE_candidates_after_MCHelper"

REPBASE_NT_DB_PREFIX="${TE_STEP_DIR}/repbase_nt"
TE_PEP_DB_PREFIX="${TE_STEP_DIR}/te_peptides"

###############################################################################
# Step 2: remove known TE consensus sequences
###############################################################################
#BLASTN thresholds used to identify TE consensus sequences already represented
# A minimum query coverage of 50% was chosen instead of the classical 80-80-80
# criterion because RepeatModeler frequently generates fragmented or partial
# consensus sequences. Since the curated library is assumed to contain more
# complete TE consensuses, a partial RepeatModeler consensus may still represent
# the same TE family even if it covers only part of the curated sequence.
# Therefore, candidates showing >=80% nucleotide identity over at least 50% of
# their own length are considered already represented in the curated library and
# are removed from further analyses.

KNOWN_TE_MIN_IDENTITY="80"
KNOWN_TE_MIN_QCOV="50"
KNOWN_TE_MIN_LENGTH="80"

###############################################################################
# Step 4: final BLAST filtering
###############################################################################
# After MCHelper, candidate sequences from all genomes are merged and
# dereplicated using CD-HIT-EST. Clustering is performed using an 80-80-80
# criterion (80% identity and at least 80% coverage of both the shorter and
# longer sequences) to collapse redundant representatives of the same TE family.
# The non-redundant candidates are then searched against RepBase (BLASTN) and a
# TE peptide database (BLASTX). Candidates with significant similarity to either
# database (E-value <= 1e-10) are retained as high-confidence TE candidates
# that will be further inspected.

# It is unlikely that bona fide canonical TEs would lack any detectable similarity 
# to previously described TE nucleotide or protein sequences. As the aim of this pipeline
# is to generate a curated library of bona fide TEs, only candidates supported by sequence 
# similarity are retained. Although some of these sequences could represent genuine
# non-canonical TEs, demonstrating their transposable nature would require extensive 
# manual investigation. 

FINAL_BLAST_EVALUE="1e-10"

CDHIT_IDENTITY="0.80"
CDHIT_SHORTER_COVERAGE="0.80"
CDHIT_LONGER_COVERAGE="0.80"

###############################################################################
# Thread settings
###############################################################################

REPEATMODELER_THREADS="64"
FILTER_THREADS="32"
MCHELPER_THREADS="64"

###############################################################################
# RepeatModeler settings
###############################################################################

# RepeatModeler option used to set the number of CPU threads.
#
# RepeatModeler versions use different option names:
#   - RepeatModeler 2.0.3 and older: use "-pa"
#   - RepeatModeler 2.0.4 and newer: use "-threads"
#
# To check your version, run:
#   RepeatModeler -h | head
#
# Then set:
#   REPEATMODELER_THREADS_OPTION="-pa"       # for RepeatModeler 2.0.3 and older
#   REPEATMODELER_THREADS_OPTION="-threads"  # for RepeatModeler 2.0.4 and newer
REPEATMODELER_THREADS_OPTION="-threads"

# Whether to run RepeatModeler with LTR structural discovery.
#
# Use "yes" to add:
#   -LTRStruct
#
# Use "no" to run RepeatModeler without LTR structural discovery.
#
# LTRStruct may improve recovery and classification of LTR retrotransposons,
# but it increases runtime and may require additional LTR-related tools in the
# RepeatModeler environment.
REPEATMODELER_USE_LTRSTRUCT="yes"

###############################################################################
# MCHelper settings
###############################################################################

# MCHelper iteratively extends each candidate sequence by a fixed number of
# nucleotides on both sides to recover potentially missing TE boundaries.
#
# MCHelper parameter "-x": Number of extension iterations.
#
# MCHelper parameter "-e": Number of nucleotides added to each side of the 
# candidate sequence at each iteration.
#
# Total maximum extension = x * e bp on each side.
#
# Default values (2 * 2000 bp) generally allow recovery of TE termini for
# fragmented RepeatModeler consensus sequences while avoiding excessive
# extension into surrounding genomic regions.

MCHELPER_X="2"
MCHELPER_EXTENSION="3000"

###############################################################################
# SLURM settings
###############################################################################
# These values are used only by slurm/run_pipeline_slurm.sh.
# Leave SLURM_PARTITION empty to use the cluster default partition.
# Leave memory or time values empty to use the cluster defaults.
#
# ARRAY_CONCURRENCY controls how many array tasks can run at the same time.
# Example: ARRAY_CONCURRENCY="3" submits arrays as 1-N%3.
# Leave it empty to let SLURM use the cluster default.

SLURM_PARTITION=""

STEP1_CPUS=""
STEP1_MEM=""
STEP1_TIME=""

STEP2_CPUS=""
STEP2_MEM=""
STEP2_TIME=""

STEP3_CPUS=""
STEP3_MEM=""
STEP3_TIME=""

STEP4_CPUS=""
STEP4_MEM=""
STEP4_TIME=""

ARRAY_CONCURRENCY=""

###############################################################################
# Software environments
###############################################################################
###############################################################################
# Software environments
###############################################################################
# Edit these functions for your cluster.
#
# Important:
# The pipeline scripts expect the commands below to be available in PATH:
#   BuildDatabase
#   RepeatModeler
#   makeblastdb
#   blastn
#   blastx
#   seqkit
#   cd-hit-est
#   python3
#
# If RepeatModeler is inside a Singularity/Apptainer container, this config
# creates temporary wrapper commands so the pipeline scripts do not need to be
# edited.

###############################################################################
# RepeatModeler environment
###############################################################################

# Path to the Dfam/TETools Singularity/Apptainer image.
TETOOLS_SIF="${SHARED:-/shared}/containers/dfam-tetools-latest.sif"

# Use "apptainer" or "singularity", depending on your cluster.
TETOOLS_CONTAINER_CMD="apptainer"

load_repeatmodeler_environment() {
    set +u

    module purge || true

    # Load the container runtime. Keep the one that exists on your cluster.
    module load apptainer || module load singularity || true

    set -u

    if ! command -v "${TETOOLS_CONTAINER_CMD}" >/dev/null 2>&1; then
        echo "ERROR: container command not found: ${TETOOLS_CONTAINER_CMD}" >&2
        echo "Edit TETOOLS_CONTAINER_CMD or load the correct module." >&2
        exit 1
    fi

    if [[ ! -s "${TETOOLS_SIF}" ]]; then
        echo "ERROR: TETools container image not found: ${TETOOLS_SIF}" >&2
        exit 1
    fi

    # Create temporary wrapper commands for the current run.
    # The pipeline scripts call BuildDatabase and RepeatModeler normally;
    # these wrappers redirect those calls into the container.
    TETOOLS_WRAPPER_DIR="${PROJECT_DIR}/.tetools_wrappers"
    mkdir -p "${TETOOLS_WRAPPER_DIR}"

    cat > "${TETOOLS_WRAPPER_DIR}/BuildDatabase" <<EOF
#!/usr/bin/env bash
exec "${TETOOLS_CONTAINER_CMD}" exec "${TETOOLS_SIF}" BuildDatabase "\$@"
EOF

    cat > "${TETOOLS_WRAPPER_DIR}/RepeatModeler" <<EOF
#!/usr/bin/env bash
exec "${TETOOLS_CONTAINER_CMD}" exec "${TETOOLS_SIF}" RepeatModeler "\$@"
EOF

    chmod +x "${TETOOLS_WRAPPER_DIR}/BuildDatabase" \
             "${TETOOLS_WRAPPER_DIR}/RepeatModeler"

    export PATH="${TETOOLS_WRAPPER_DIR}:${PATH}"
}

###############################################################################
# Filtering environment
###############################################################################

load_filtering_environment() {
    set +u

    module purge || true

    if [[ -n "$(_container_images "$ENV_TE_FILTER")" ]]; then
        # $CONTAINER_DIR/te_filter/ (BLAST+, cd-hit, seqkit images - containers/pull_containers.sh)
        activate_env "$ENV_TE_FILTER"
    else
        module load ncbi-blast/2.12.0
        module load seqkit/2.9.0

        # Keep the module name that exists on your cluster.
        module load cd-hit || module load cd-hit-est || true
    fi

    set -u

    echo "Filtering environment:"

    if ! command -v makeblastdb >/dev/null 2>&1; then
        echo "ERROR: makeblastdb not found after loading filtering modules." >&2
        exit 1
    fi

    if ! command -v blastn >/dev/null 2>&1; then
        echo "ERROR: blastn not found after loading filtering modules." >&2
        exit 1
    fi

    if ! command -v blastx >/dev/null 2>&1; then
        echo "ERROR: blastx not found after loading filtering modules." >&2
        exit 1
    fi

    if ! command -v seqkit >/dev/null 2>&1; then
        echo "ERROR: seqkit not found after loading filtering modules." >&2
        exit 1
    fi

    if ! command -v cd-hit-est >/dev/null 2>&1; then
        echo "ERROR: cd-hit-est not found after loading filtering modules." >&2
        exit 1
    fi

    echo "  makeblastdb: $(command -v makeblastdb)"
    echo "  blastn:      $(command -v blastn)"
    echo "  blastx:      $(command -v blastx)"
    echo "  seqkit:      $(command -v seqkit)"
    echo "  cd-hit-est:  $(command -v cd-hit-est)"
}

###############################################################################
# MCHelper environment
###############################################################################

load_mchelper_environment() {
    set +u

    module purge || true

    activate_env "$ENV_MCHELPER"  # env names / backend: pipeline_config.sh
    set -u

    if ! command -v mchelper >/dev/null 2>&1 && ! command -v python3 >/dev/null 2>&1; then
        echo "ERROR: neither mchelper nor python3 found after activating MCHelper environment." >&2
        exit 1
    fi
}
