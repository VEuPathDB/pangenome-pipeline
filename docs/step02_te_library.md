# Step 02 – TE-library pipeline (SLURM)

[`workflow/step02_te_annotation_clean_up/te_library/`](../workflow/step02_te_annotation_clean_up/te_library/)
holds an independent SLURM version of the TE candidate discovery pipeline described in
[step 02](step02_te_annotation_clean_up.md#2---te-candidate-discovery-pipeline).

These SLURM scripts were adapted for the cluster environment used during testing.
Software modules, micromamba paths, container paths, environment names, file
paths, and SLURM resource settings may not match other HPC systems.

Before running the pipeline on another cluster, users should review and edit the
files in `te_library/` as needed.

## Files

```text
workflow/step02_te_annotation_clean_up/te_library/
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
$PIPELINE_DIR/workflow/step02_te_annotation_clean_up/te_library/config.sh
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

The cluster paths themselves (MCHelper, RepBase, TE peptides, BUSCO lineage, micromamba) are set
in `pipeline_config.sh` — see the [main README](../README.md#getting-the-pipeline-into-a-project).

Edit SLURM resources in the header of each `.slurm.sh` script:

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
```

## Where things go

Every script is run **from the project folder** (the one with `genomes/`), like every other job:

- inputs: `genomes/*.fna`, and the known-TE library `CURATED_TE_LIBRARY` (default
  `$STEP02_OUT/fake_lib.fasta`);
- outputs: `$STEP02_OUT/te_library/` — `config/` (the genome and MCHelper tables),
  `repeatmodeler_more_genomes/`, `TE_candidate_screen_original_headers/`, `MCHelper_newTEs/`,
  `final_TE_candidates_after_MCHelper/` and the BLAST databases;
- logs: `logs/` in the project folder (the runner's own log too). Create it once with
  `mkdir -p logs` — SLURM does not create it, and a job without it leaves no log.

`$STEP02_OUT` is `$OUTPUT_DIR/step02_te_annotation_clean_up`, by default
`output/step02_te_annotation_clean_up/` in the project folder (see
[outputs](../README.md#outputs)).

## How to Run

**First, after copying the pipeline to the cluster**, make the scripts executable — Google
Drive and some transfer tools drop the executable bit, which makes Step 2 fail with
`Permission denied`:

```bash
chmod +x $PIPELINE_DIR/workflow/step02_te_annotation_clean_up/te_library/*.sh
```

The helpers run with `bash` find `config.sh` next to themselves. The `StepN_*.slurm.sh` jobs are
copied by `sbatch`, so they find the pipeline through `PIPELINE_DIR` instead (set in your
`~/.bashrc`, see the main README; `run_pipeline_slurm.sh` passes it on to the jobs it submits).

### Option 1: run all steps sequentially

```bash
cd <project>
bash $PIPELINE_DIR/workflow/step02_te_annotation_clean_up/te_library/run_pipeline_slurm.sh
```

This runner uses `sbatch --wait`. It submits one SLURM job, waits until that job
finishes, and then starts the next step. If a step fails, the runner stops.
It runs `check_dependencies.sh` before preparing the genome table.

### Option 2: run each step manually

From the project folder:

```bash
TE=$PIPELINE_DIR/workflow/step02_te_annotation_clean_up/te_library

bash $TE/check_dependencies.sh
bash $TE/00_prepare_genome_table.sh

sbatch $TE/Step1_run_repeatmodeler.slurm.sh
# wait for it to finish:
#   squeue -u "$USER"
#   tail -n 80 logs/TE_step1_repeatmodeler_*.{out,err}

sbatch $TE/Step2_run_TE_candidates_original_headers.slurm.sh
bash $TE/02_prepare_mchelper_table.sh

sbatch $TE/Step3_run_mchelper.slurm.sh
sbatch $TE/Step4_run_final_TE_filter_after_MCHelper.slurm.sh
```

## Resume Behaviour

The scripts skip outputs that already exist (all under `$STEP02_OUT/te_library/`):

- Step 1 skips a genome if `consensi.fa.classified` already exists.
- Step 2 skips if `TE_candidate_screen_original_headers/summary.tsv` exists.
- Step 3 skips a genome if `curated_sequences_NR.fa` already exists.
- Step 4 skips if `final_TE_candidates_after_MCHelper/final_potential_new_TEs.fa`
  exists.

To rerun a step, remove the corresponding output folder or file before
submitting the step again.

[Back to step 02](step02_te_annotation_clean_up.md)
