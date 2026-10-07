# Running the Pipeline with SLURM

This folder contains an independent SLURM version of the TE candidate discovery
pipeline.

These SLURM scripts were adapted for the cluster environment used during testing.
Software modules, micromamba paths, container paths, environment names, file
paths, and SLURM resource settings may not match other HPC systems.

Before running the pipeline on another cluster, users should review and edit the
files in the `slurm/` folder as needed.

If you run the SLURM version, edit and use only the files inside `slurm/`.
The non-SLURM scripts in the main folder are not required.

## Files

```text
slurm/
├── README.md
├── config.sh
├── check_dependencies.sh
├── 00_prepare_genome_table.sh
├── Step1_run_repeatmodeler.slurm.sh
├── Step2_run_TE_candidates_original_headers.slurm.sh
├── TE_candidate_pipeline_original_headers.sh
├── 02_prepare_mchelper_table.sh
├── Step3_run_mchelper.slurm.sh
├── Step4_run_final_TE_filter_after_MCHelper.slurm.sh
└── run_pipeline_slurm.sh
```

## What to Edit

Edit the SLURM configuration file:

```bash
02_TE_annotation_Clean_Up/slurm/config.sh
```

This file contains:

- input genome pattern;
- curated TE library path;
- BUSCO lineage path;
- MCHelper path;
- RepBase nucleotide FASTA path;
- TE peptide FASTA path;
- RepeatModeler, filtering, MCHelper, and CD-HIT parameters;
- cluster-specific software loading commands.

Edit SLURM resources in the header of each `.slurm.sh` script:

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
```

## How to Run

**First, after copying the folder to the cluster**, make the scripts executable — Google
Drive and some transfer tools drop the executable bit, which makes Step 2 fail with
`Permission denied`:

```bash
chmod +x 02_TE_annotation_Clean_Up/slurm/*.sh
```

Two ways the scripts find the pipeline folder — they are **not** the same:

- `run_pipeline_slurm.sh`, `check_dependencies.sh`, `00_prepare_genome_table.sh`,
  `02_prepare_mchelper_table.sh` locate themselves and can be run **from any directory**.
- The `StepN_*.slurm.sh` jobs are copied by `sbatch`, so they **cannot** locate themselves.
  They resolve `slurm/config.sh` and their `logs/` folder from the directory you run
  `sbatch` in. Submit them **from inside `02_TE_annotation_Clean_Up/`**, or they fail with
  `config.sh: No such file`.

`run_pipeline_slurm.sh` handles this for you (it `cd`s into `02_TE_annotation_Clean_Up/`
before every `sbatch`). It only matters when you re-run a single step by hand.

### Option 1: run all steps sequentially

```bash
bash 02_TE_annotation_Clean_Up/slurm/run_pipeline_slurm.sh
```

This runner uses `sbatch --wait`. It submits one SLURM job, waits until that job
finishes, and then starts the next step. If a step fails, the runner stops.
It runs `check_dependencies.sh` before preparing the genome table.

### Option 2: run each step manually

The self-locating helpers can be run from the project root:

```bash
bash 02_TE_annotation_Clean_Up/slurm/check_dependencies.sh
bash 02_TE_annotation_Clean_Up/slurm/00_prepare_genome_table.sh
```

The `sbatch` steps must be submitted **from inside `02_TE_annotation_Clean_Up/`**:

```bash
cd 02_TE_annotation_Clean_Up

sbatch slurm/Step1_run_repeatmodeler.slurm.sh
# wait for it to finish:
#   squeue -u "$USER"
#   tail -n 80 logs/slurm/TE_step1_repeatmodeler_*.{out,err}

sbatch slurm/Step2_run_TE_candidates_original_headers.slurm.sh
# then, back from the project root or from here:
bash slurm/02_prepare_mchelper_table.sh

sbatch slurm/Step3_run_mchelper.slurm.sh
sbatch slurm/Step4_run_final_TE_filter_after_MCHelper.slurm.sh
```

## Resume Behaviour

The scripts skip outputs that already exist:

- Step 1 skips a genome if `consensi.fa.classified` already exists.
- Step 2 skips if `TE_candidate_screen_original_headers/summary.tsv` exists.
- Step 3 skips a genome if `curated_sequences_NR.fa` already exists.
- Step 4 skips if `final_TE_candidates_after_MCHelper/final_potential_new_TEs.fa`
  exists.

To rerun a step, remove the corresponding output folder or file before
submitting the step again.
