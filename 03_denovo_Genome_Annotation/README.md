# 03 – *de novo* Genome Annotation with Funannotate

> ⚠️ **Not reviewed yet** — draft content, adapted from Eve's pipeline and my notes

***Purpose:*** Comprehensive guide and SLURM pipeline for annotating fungal genomes *de novo* using **Funannotate** via the shared Singularity installation on the HPC cluster.

Skip this step if good gene annotations are already available for all your genomes. 

---

## Table of Contents

- [Overview](#overview)
- [HPC Environment Setup & Prerequisites](#hpc-environment-setup--prerequisites)
  - [1. Environment variables (`~/.bashrc`)](#1-environment-variables-bashrc)
  - [2. GeneMark academic license key](#2-genemark-academic-license-key)
  - [3. PASA configuration](#3-pasa-configuration)
  - [4. Verifying dependencies and installation](#4-verifying-dependencies-and-installation)
- [Step 1: Genome Assembly Preparation](#step-1-genome-assembly-preparation)
  - [Clean & Sort contigs](#clean--sort-contigs)
  - [Repeat Masking](#repeat-masking)
- [Step 2: RNA-seq Evidence Preparation](#step-2-rna-seq-evidence-preparation)
  - [Download RNA-seq data from SRA/ENA](#download-rna-seq-data-from-sraena)
  - [Quality Control](#quality-control)
  - [Combining reads & fixing FASTQ headers](#combining-reads--fixing-fastq-headers)
  - [Read Normalization for SRA Paired-End Data](#read-normalization-for-sra-paired-end-data)
- [Step 3: Funannotate Train](#step-3-funannotate-train)
- [Step 4: Funannotate Predict](#step-4-funannotate-predict)
- [Step 5: Funannotate Update (PASA UTRs & Alternative Splicing)](#step-5-funannotate-update-pasa-utrs--alternative-splicing)
- [Step 6: Funannotate Annotate (Functional Annotation)](#step-6-funannotate-annotate-functional-annotation)
- [Step 7: Quality Assessment & BUSCO](#step-7-quality-assessment--busco)
- [Troubleshooting & Practical Notes](#troubleshooting--practical-notes)

---

## Overview

[Funannotate](https://funannotate.readthedocs.io/) is an automated genome annotation pipeline designed for fungi and other small eukaryotes. The pipeline integrates evidence across multiple stages:

```
[ Genome Assembly ] ──► funannotate clean & sort ──► RepeatMasker (Soft-masking)
                                                            │
[ RNA-seq Reads   ] ──► FastQC ──► fix_headers.py           │
                                       │                    ▼
                                       └──────────► funannotate train (Trinity + PASA)
                                                            │
                                                            ▼
                                                    funannotate predict (EVM + Ab initio)
                                                            │
                                                            ▼
                                                    funannotate update (PASA / UTRs)
                                                            │
                                                            ▼
                                                    funannotate annotate (Pfam, GO, InterProScan)
                                                            │
                                                            ▼
                                                    BUSCO Evaluation & Final GFF3/GBK
```

**Shared Cluster Resources:**
- **Shared Directory**: `/mnt/hc-storage/groups/veupathdb_annotation/eb_tools/funannotate`
- **Main Container**: `${SHARED}/containers/funannotate.sif` (v1.8.17)
- **Repeat Masking Container**: `${SHARED}/containers/dfam-tetools-latest.sif`
- **Databases**: `${SHARED}/databases/funannotate_db/` (UniProt, Pfam, MEROPS, dbCAN, GO, BUSCO)

[↑ Back to top](#table-of-contents)

---

## HPC Environment Setup & Prerequisites

### 1. Environment variables (`~/.bashrc`)

Set `FUNANNOTATE_SHARED` (the shared funannotate install) in the project's `pipeline_config.sh`
(see the [main README](../README.md#getting-the-pipeline-into-a-project)) - the job scripts read it
from there. For running funannotate by hand, and for your personal Singularity cache directory, add
these to your `~/.bashrc`:

```bash
# Funannotate shared installation
export FUNANNOTATE_SHARED=/mnt/hc-storage/groups/veupathdb_annotation/eb_tools/funannotate
export FUNANNOTATE_DB=${FUNANNOTATE_SHARED}/databases/funannotate_db

# Singularity / Apptainer cache — MUST point to your own user space (never /tmp)
export SINGULARITY_CACHEDIR=/mnt/hc-storage/users/${USER}/singularity_cache
export APPTAINER_CACHEDIR=/mnt/hc-storage/users/${USER}/singularity_cache
export SINGULARITY_TMPDIR=/mnt/hc-storage/users/${USER}/singularity_cache
export APPTAINER_TMPDIR=/mnt/hc-storage/users/${USER}/singularity_cache
```

Apply changes and create cache:
```bash
source ~/.bashrc
mkdir -p "${SINGULARITY_CACHEDIR}"
```

### 2. GeneMark academic license key

GeneMark-ES is included in the shared installation, but each user must obtain their own free academic license key:

1. Visit the [GeneMark License Download](http://exon.gatech.edu/GeneMark/license_download.cgi) page.
2. Select **GeneMark-ES/ET/EP+ ver 4.73_lic for LINUX 64 kernel 3.10-5**.
3. Download the key (`gm_key_64.gz`), copy it to your HPC home directory, and uncompress:
   ```bash
   gunzip ~/gm_key_64.gz
   cp ~/gm_key_64 ~/.gm_key
   ```

### 3. PASA configuration

The shared installation uses an SQLite-based PASA configuration (`${FUNANNOTATE_SHARED}/pasa.config.txt`). When running with SQLite, **no manual database configuration or user-specific editing of `pasa.config.txt` is required**.

### 4. Verifying dependencies and installation

Test that the container and dependencies are functioning:

```bash
singularity exec \
    -B "${FUNANNOTATE_SHARED}/databases/funannotate_db:${FUNANNOTATE_SHARED}/databases/funannotate_db" \
    -B "${FUNANNOTATE_SHARED}/overlay/usr/local/bin:/usr/local/bin" \
    -B "${FUNANNOTATE_SHARED}/overlay/usr/lib/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu" \
    --env FUNANNOTATE_DB="${FUNANNOTATE_SHARED}/databases/funannotate_db" \
    --env GENEMARK_PATH="/usr/local/bin/genemark" \
    --env PATH="/usr/local/bin:/venv/bin:/usr/bin:/bin" \
    "${FUNANNOTATE_SHARED}/containers/funannotate.sif" \
    funannotate check --show-versions
```

[↑ Back to top](#table-of-contents)

---

## Step 1: Genome Assembly Preparation

### Clean & Sort contigs

Removes small duplicate contigs (<500 bp) and standardizes FASTA headers.

Script: [`scripts/clean_and_sort.slurm.sh`](scripts/clean_and_sort.slurm.sh)

```bash
sbatch scripts/clean_and_sort.slurm.sh
```

### Repeat Masking

Soft-mask repetitive elements using RepeatMasker prior to gene prediction. Soft-masking converts repetitive sequences to lowercase nucleotides (allowing aligners and gene predictors to utilize these regions when supported by strong transcript/protein evidence, while suppressing false-positive *ab initio* gene calls inside TEs).

#### Choosing the Repeat Library:
- **Curated TE Library from Step 02 (Recommended):** If **[02 – TE Annotation and Clean Up](../02_TE_annotation_Clean_Up/README.md)** has already been conducted for the target species/genus, use the final **curated TE library** (`curated_TE_library.fa`) produced in that step. Using a manually curated library avoids over-masking conserved host genes and minimizes annotation artifacts.
- **RepeatModeler *de novo* Library (Alternative):** If Step 02 has not yet been performed or for standalone annotations, you can use the consensus output library generated *de novo* by RepeatModeler (e.g., `*-families.fa`) or a relevant clade library from Repbase/Dfam.

Script: [`scripts/repeatmask.slurm.sh`](scripts/repeatmask.slurm.sh)

```bash
sbatch scripts/repeatmask.slurm.sh
```

[↑ Back to top](#table-of-contents)

---

## Step 2: RNA-seq Evidence Preparation

### Download RNA-seq data from SRA/ENA

Download selected RNA-seq accessions (SRR/ERR) using `prefetch` and convert to FASTQ with `fasterq-dump`:

```bash
module load sra-tools/3.0.3

# Download .sra files
prefetch --option-file accessions.txt -O sra_files --max-size 50G

# Convert to FASTQ
mkdir -p fastq
while read -r acc; do
    fasterq-dump --split-files --threads 8 -O fastq "sra_files/${acc}/${acc}.sra"
done < accessions.txt

gzip fastq/*.fastq
```

Verify paired-end integrity (both `_1.fastq.gz` and `_2.fastq.gz` must have identical read counts):
```bash
seqkit stats -T -j 4 fastq/*.fastq.gz
```

### Quality Control

Script: [`scripts/fastqc.slurm.sh`](scripts/fastqc.slurm.sh)

```bash
sbatch scripts/fastqc.slurm.sh
```

### Combining reads & fixing FASTQ headers

1. **Concatenate all runs:**
   ```bash
   cat fastq/*_1.fastq.gz > combined_left_trinity.fastq.gz
   cat fastq/*_2.fastq.gz > combined_right_trinity.fastq.gz
   ```

2. **Fix FASTQ headers for Trinity:**
   Trinity requires standard `/1` and `/2` read identifiers and strict `+` delimiters on the third line.

   Scripts: [`scripts/fix_headers.py`](scripts/fix_headers.py) & [`scripts/fix_headers.slurm.sh`](scripts/fix_headers.slurm.sh)
   ```bash
   sbatch scripts/fix_headers.slurm.sh
   ```

### Read Normalization for SRA Paired-End Data

> [!WARNING]
> **Trinity Normalization Bug with SRA Reads:**
> By default, `funannotate train` internally triggers Trinity's `insilico_read_normalization.pl`. However, when using paired-end RNA-seq reads downloaded from SRA/ENA, Trinity's internal Perl merge script (`nbkc_merge_left_right_stats.pl`) fails to match read pairs due to a lexicographical sorting mismatch caused by the `_1/1` vs `_2/2` suffixes.

To handle this, use one of the two strategies:

1. **Option A — Normalize with `merge_paired_stats.py` (Standard for pooled/large datasets):**
   - If Trinity normalization fails on k-mer merging, navigate to the `tmp_normalized_reads/` directory.
   - Run [`scripts/merge_paired_stats.py`](scripts/merge_paired_stats.py) to merge `left.fa.K25.stats.sort` and `right.fa.K25.stats.sort` into `pairs.K25.stats.custom`.
   - Finish extracting `left.norm.fq` and `right.norm.fq`.
   - Provide the normalized reads directly to `funannotate train` with `--left left.norm.fq --right right.norm.fq --no_normalize_reads`.

2. **Option B — Skip normalization:**
   - If your RNA-seq dataset is small (<20–30 million reads) or has already been normalized, pass `--no_normalize_reads` directly to `funannotate train`.

[↑ Back to top](#table-of-contents)

---

## Step 3: Funannotate Train

Executes Trinity *de novo* transcript assembly, aligns transcripts back to the genome with GMAP, and constructs training gene models with PASA.

Script: [`scripts/train.slurm.sh`](scripts/train.slurm.sh)

```bash
sbatch scripts/train.slurm.sh
```

**Key Outputs in `funannotate_train/`:**
- `funannotate_train/training/trinity.fasta.clean`: Cleaned Trinity transcript assemblies.
- `funannotate_train/training/pasa/`: PASA alignment and transcript assemblies.
- `funannotate_train/training/left.fq.gz` & `right.fq.gz`: Cleaned/normalized reads passed to downstream steps.

[↑ Back to top](#table-of-contents)

---

## Step 4: Funannotate Predict

Combines evidence from `funannotate train`, protein alignments (Exonerate), and *ab initio* gene predictors (GeneMark-ES, AUGUSTUS, SNAP, GlimmerHMM) using EvidenceModeler (EVM).

Script: [`scripts/predict.slurm.sh`](scripts/predict.slurm.sh)

```bash
sbatch scripts/predict.slurm.sh
```

> **Augustus Configuration Note:** Augustus requires a writable config directory. `predict.slurm.sh` automatically copies the Augustus configuration into your local working directory (`augustus_config_predict/`) and bind-mounts it.

**Key Outputs in `funannotate_predict/predict_results/`:**
- `<Species>_<Strain>.gff3`: Raw predicted gene models.
- `<Species>_<Strain>.proteins.fa`: Translated protein sequences.
- `<Species>_<Strain>.gbk`: GenBank file of predictions.

[↑ Back to top](#table-of-contents)

---

## Step 5: Funannotate Update (PASA UTRs & Alternative Splicing)

Passes the EVM predictions through PASA with the RNA-seq data to add 5' and 3' UTRs, refine exon junctions, and detect alternative splicing isoforms.

Script: [`scripts/update.slurm.sh`](scripts/update.slurm.sh)

```bash
sbatch scripts/update.slurm.sh
```

**Key Outputs in `funannotate_update/update_results/`:**
- `<Species>_<Strain>.gff3`: Updated gene models with annotated UTRs and isoforms.
- `<Species>_<Strain>.proteins.fa`: Protein sequences from refined models.
- `<Species>_<Strain>.tbl`: NCBI Sequin table file.

[↑ Back to top](#table-of-contents)

---

## Step 6: Funannotate Annotate (Functional Annotation)

Performs comprehensive functional annotation, assigning Pfam domains, GO terms, UniProt/SwissProt IDs, MEROPS proteases, dbCAN CAZymes, SignalP peptides, and InterProScan domains.

Script: [`scripts/annotate.slurm.sh`](scripts/annotate.slurm.sh)

```bash
sbatch scripts/annotate.slurm.sh
```

*Optional standalone InterProScan:* If running InterProScan independently outside the container:
Script: [`scripts/iprscan.slurm.sh`](scripts/iprscan.slurm.sh)

```bash
sbatch scripts/iprscan.slurm.sh
```

[↑ Back to top](#table-of-contents)

---

## Step 7: Quality Assessment & BUSCO

Evaluates annotation completeness using BUSCO against appropriate lineage datasets (e.g. `fungi_odb10`, `ascomycota_odb10`, or `dikarya`).

Script: [`scripts/busco.slurm.sh`](scripts/busco.slurm.sh)

```bash
sbatch scripts/busco.slurm.sh
```

[↑ Back to top](#table-of-contents)

---

## Troubleshooting & Practical Notes

Based on practical testing and optimization notes (`Notes_AL.Rmd`):

1. **Singularity PATH Order for Exonerate:**
   The container's internal `/venv/bin` may contain an unpatched Exonerate binary. Ensure `/usr/local/bin` precedes `/venv/bin` in container `PATH`:
   ```bash
   --env PATH="/usr/local/bin:/venv/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
   ```

2. **GeneMark Path:**
   Always pass `--env GENEMARK_PATH="/usr/local/bin/genemark"` in the Singularity command.

3. **Augustus Config Permissions:**
   Augustus fails if it cannot write training parameters to its config directory. Always copy `/usr/share/augustus/config` to your local working directory and mount it:
   ```bash
   -B "${PWD}/augustus_config_predict:/usr/share/augustus/config"
   ```

4. **Java Heap Memory:**
   PASA and Trinity steps can require substantial Java memory. Set `export _JAVA_OPTIONS="-Xmx32g"` prior to executing `funannotate update` or `funannotate train`.

5. **Trinity FASTQ Third Line:**
   Trinity throws errors if the 3rd line of FASTQ entries contains redundant header text instead of just `+`. Always preprocess pooled FASTQ files using [`scripts/fix_headers.py`](scripts/fix_headers.py).

[↑ Back to top](#table-of-contents) &nbsp;|&nbsp; [Back to top level](../README.md)
