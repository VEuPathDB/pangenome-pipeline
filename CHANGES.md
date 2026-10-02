# Changes on the `add-python-tooling` branch

Everything that changed compared with `main`: 49 commits, made in this order:

1. Python tooling and tests
2. Removing duplicated code
3. Removing hardcoded paths
4. Moving the Python into a package
5. The new `workflow/` + `docs/` + `output/` layout
6. Running without micromamba

Each commit message has the full detail (`git log main..add-python-tooling`). This page is the
overview for anyone testing the branch.

**Contents:** [Quick start](#quick-start-for-testing) ·
[Layout](#1-layout-where-things-are-now) · [Outputs](#2-outputs) ·
[Configuration](#3-configuration-instead-of-hardcoded-values) ·
[Tool environments](#4-tool-environments-micromamba-is-now-optional) ·
[Python code](#5-python-code) ·
[Behaviour changes](#6-changes-that-can-change-results) ·
[Tooling and tests](#7-tooling-and-tests) · [Not yet verified](#8-not-yet-verified--open-items)

---

## Quick start for testing

```bash
# 1. The pipeline is its own checkout, separate from the species project
git clone <repo> pangenome-pipeline && cd pangenome-pipeline
git switch add-python-tooling
export PIPELINE_DIR=$PWD

# 2. Python for the pipeline's own scripts, pinned by uv.lock (no conda needed)
uv sync --frozen --no-dev

# 3. Cluster settings
cp pipeline_config.example.sh pipeline_config.sh     # fill in what you use

# 4. Optional: tool containers instead of micromamba (see section 4)
bash containers/pull_containers.sh /shared/path/containers   # then set CONTAINER_DIR

# 5. Jobs are submitted from the species project folder
cd /path/to/my_species_project          # has genomes/, genome_order.txt, ...
sbatch $PIPELINE_DIR/workflow/step01_genome_selection_qc/synteny_vs_ref.slurm.sh
```

Projects started before this branch: follow [docs/migrating_projects.md](docs/migrating_projects.md),
which has a script that moves existing results to their new places.

---

## 1. Layout: where things are now

| What | Where |
|---|---|
| SLURM job scripts | `workflow/stepNN_<name>/` (e.g. `workflow/step05_cross_genome_annotation_integration/`) |
| Step documentation | `docs/stepNN_<name>.md` |
| Python code | `src/pangenome_pipeline/stepNN_<name>/`, shared helpers in `src/pangenome_pipeline/common/` |
| Shared shell helpers | `lib/slurm_common.sh` (every job loads it) |
| Species-only scripts | `species_specific/<Species>/stepNN_<name>/` |
| Get-PanLoci alternative | `workflow/alternative_get_panloci/` (+ `no_slurm/`), `docs/alternative_get_panloci.md` |
| Container recipes | `containers/` |
| Species reports | `results/` (unchanged) |

- The step numbers are kept. Python package names can't start with a digit, so folders are
  called `step01_...`, `step02_...` and so on.
- The pipeline checkout and the species project are now **separate folders**. Jobs are submitted
  from the project folder and find the pipeline through `$PIPELINE_DIR`, or from where the script
  sits.
- Python runs as `python3 -m pangenome_pipeline.<step>.<module>`. `lib/slurm_common.sh` puts
  `src/` on `PYTHONPATH`, so nothing has to be installed.
- The full old → new tables (pipeline files, results, logs) are in
  [docs/migrating_projects.md](docs/migrating_projects.md).

## 2. Outputs

- Every step writes under `$OUTPUT_DIR/stepNN_<name>/`. `OUTPUT_DIR` defaults to `output/` in
  the project folder.
- In job scripts these folders are `$STEP01_OUT` … `$STEP08_OUT` and `$ALT_PANLOCI_OUT`. Each can
  be pointed elsewhere on its own.
- `OUTPUT_DIR` can be changed in `pipeline_config.sh` or per job, e.g.
  `sbatch --export=ALL,OUTPUT_DIR=/scratch/run2 ...`. A relative path is taken from the project
  folder.
- SLURM logs go to `logs/<job>_<id>.out` / `.err` in the project.
- Step 03 runs inside `$STEP03_OUT`, which takes the role of its old working folder (inputs
  copied in, all funannotate output there).
- Get-PanLoci: `RUN_DIR` defaults to `$ALT_PANLOCI_OUT` instead of a placeholder path that had
  to be edited.
- Step 06: get_pangenes is run from its input folder's parent, so its `<input>_pangenes/` output
  lands in `$STEP06_OUT`, not in the submit folder.

## 3. Configuration instead of hardcoded values

**`pipeline_config.sh`.** All cluster-specific paths moved here, out of the scripts.
- Start from `pipeline_config.example.sh`. The file isn't committed.
- The pipeline's copy is loaded first, then the project's own `pipeline_config.sh` if it has one,
  so a project only sets what differs.
- Any setting can also be given per job with `sbatch --export=ALL,NAME=value`.
- A job that needs an empty setting stops with a message naming it.
- Settings: `OUTPUT_DIR`, `ENV_BACKEND`, `CONTAINER_DIR`, `CONTAINER_BINDS`, `PIPELINE_VENV`,
  `MAMBA_EXE` / `MAMBA_ROOT_PREFIX`, the environment names (`ENV_SYNTENY`, `ENV_BIOINFO`, …),
  `TE_AID`, `MCHELPER_PY`, `REPBASE_NT_FASTA`, `TE_PEP_FASTA`, `BUSCO_LINEAGE`,
  `FUNANNOTATE_SHARED`, `GET_PANGENES_ENV`, `INTERPRO_ENTRY_LIST`.

**Species defaults removed.** These were silently set to one species' values (e.g. *A. fumigatus*
Af293, Toxo7b). They now have to be set, and the job stops if they aren't:

| Setting | Used by | Was |
|---|---|---|
| `SPECIES`, `STRAIN` | step 03 annotate / predict / train | A. fumigatus Af293 |
| `SOURCE_GENOME`, `ORGANELLE_CONTIGS` | step 01 organelle screen | AF293's mitochondrion |
| `LABEL` | steps 07 and 08 jobs | `Toxo7b` |
| `ID_PREFIX` | step 05 `tag_final_annotation_job.sh` | `Afu` |

**Other hardcoding removed:**
- Personal and cluster paths are gone from the code. `tests/test_no_hardcoded_paths.py` keeps
  them out.
- `module load python/3.11.9` in three jobs is replaced by `activate_env`.
- The `#SBATCH --partition` / account lines can be overridden for every job without editing
  scripts: `export SBATCH_PARTITION=... SBATCH_ACCOUNT=...` before `sbatch`.

## 4. Tool environments: micromamba is now optional

Jobs call `activate_env <name>`. With `ENV_BACKEND=auto` (the default) it uses, **per
environment**, the first of these that is set up:

1. **container**: `$CONTAINER_DIR/<name>.sif`, or a folder `$CONTAINER_DIR/<name>/` of images.
   - Each program on the image's PATH gets a small wrapper that runs it with apptainer or
     singularity.
   - The project, the pipeline, `$OUTPUT_DIR` and `CONTAINER_BINDS` are visible inside.
   - Programs the host already has in `/usr/bin` (awk, perl, python3, …) are left to the host.
2. **micromamba**: `MAMBA_EXE` and `MAMBA_ROOT_PREFIX` set (the old behaviour).
3. **conda**: `conda` on PATH.
4. **path**: whatever is already loaded (`module load`, `~/.bashrc`).

The pipeline's own Python (pandas, matplotlib, numpy, scipy, seaborn) runs from the **uv venv**
pinned by `uv.lock`. Its `python3` is used whenever the environment doesn't bring one.

Apptainer (formerly Singularity) is the HPC container runtime. It needs no daemon and no root,
runs inside a SLURM job like any other program, and keeps the job's CPU, memory and time limits.

**`containers/pull_containers.sh [DIR]`** fetches pinned images into that layout. It skips ones
already present and keeps its temporary space inside `DIR`.

| Environment | Image(s) | Used by |
|---|---|---|
| `synteny/` | biocontainers minimap2, skani, samtools, blast | step 01 QC, step 02/05/08 Python jobs |
| `bioinfo/` | biocontainers agat, samtools, seqkit | steps 02, 04, 05 |
| `repeatmasker4.2.sif` | `dfam/tetools:2.00` (RepeatMasker 4.2.4 + Dfam 4.0) | step 02 RepeatMasker |
| `te_annot.sif` | biocontainers te-aid (with R + BLAST+) | step 02 TE-Aid |
| `get_pangenes.sif` | biocontainers get_pangenes | step 06 |
| `te_filter/` | biocontainers blast, cd-hit, seqkit | step 02 TE library Steps 2 and 4 |
| `mchelper.sif` | built from `containers/mchelper.def` | step 02 TE library Step 3 |

- **`containers/mchelper.def`**: MCHelper has no published image. This recipe follows MCHelper's
  own Dockerfile, with its databases and Pfam35 set up at build time. Build it with
  `APPTAINER_TMPDIR=/big/disk apptainer build --fakeroot mchelper.sif containers/mchelper.def`.
- **Using the old setup still works.** Jobs that gained a container route keep their old
  behaviour when the old setting is filled in:
  - step 06 uses `GET_PANGENES_ENV` if set, else `activate_env get_pangenes`
  - TE-Aid uses `TE_AID` if set, else `TE-Aid` from the environment
  - MCHelper uses `MCHELPER_PY` if set, else the image's `mchelper`
  - TE library filtering uses `$CONTAINER_DIR/te_filter/` if present, else its `module load`s
- **Not covered:** step 03 (funannotate), InterProScan and BUSCO still use the cluster's shared
  containers and modules.

## 5. Python code

- All Python is in the `src/pangenome_pipeline` package. No logic was rewritten in the move.
- **Duplicates removed:**
  - byte-identical copies of the Get-PanLoci scripts
  - the four step 06 wrappers now share `get_pangenes_common.sh` (wrappers kept as separate files)
- **Shared helpers, instead of a copy per script:**
  - `common/gff.py`: GFF finding, attribute parsing, biotypes, descriptions, and the
    `ID_RE` / `PARENT_RE` / `DBXREF_RE` / … patterns (were copied into 12 modules)
  - `common/intervals.py`: interval merging and overlap
  - `common/io.py`: genome lists, ID lists, PAF reading, `load_extra_genes`
  - `step05_cross_genome_annotation_integration/annotation_common.py`: step 05's parsing and overlap logic
  - `alternative_get_panloci/panloci_common.py`: the Get-PanLoci helpers
- **Deliberately left separate**, because the "duplicates" behave differently and merging them
  would change outputs: `chrom_label`, `read_fai`, `gff_seqids`, step 05's two `parse_gff`.
- Formatted and linted (ruff, isort, mypy). Invalid `'\>'` string escapes in Get-PanLoci step 1
  were fixed (same printed text).

## 6. Changes that can change results

These are bug fixes, but outputs can differ from a run on `main`:

| Step | Change |
|---|---|
| 01 | `contig_redundancy.slurm.sh` never passed `--ref`, so for any reference other than AF293 every contig showed 0 % placement with no error. It now uses the first genome in `genome_order.txt`. |
| 02 | Gene descriptions now decode `+` as a space (VEuPathDB's `hypothetical+protein` → `hypothetical protein`), as steps 01 and 05 already did. |
| 01, 02 | The description is the first *non-empty* description attribute. Trailing whitespace / `\r` is stripped from attribute values. |
| 05 | Gene renumbering broke ties (genes starting at the same position, common for Liftoff extra copies) in hash order, so IDs could swap between runs. Ties are now broken by the original gene ID. |
| 05 | `ID=`, `Parent=`, `extra_copy_number=`, `valid_ORF=` are matched as whole attributes only; e.g. `sequence_ID=` was being read as the ID if it came first. No change on Liftoff's own output order. |
| 05 | Liftoff `coverage=1e-05` was read as `1.0` and passed the `--min-rescued-coverage` floor. Coverage is now parsed in full. |

## 7. Tooling and tests

- `pyproject.toml` declares the dependencies. `uv.lock` pins them, and a test fails if the lock
  goes stale.
- `noxfile.py`: `nox -s tests` runs ruff, format, isort, mypy, toml checks and pytest.
- **312 tests**, including:
  - smoke tests: every script compiles and parses its CLI, and every shell script has valid syntax
  - unit tests for the shared helpers and step code
  - step 05 end-to-end snapshot tests (`UPDATE_SNAPSHOTS=1` regenerates them)
  - job scripts run from a fake project folder
- **Guard tests**, so the clean-up can't quietly come undone:
  - no hardcoded paths
  - no references to scripts that don't exist
  - no old folder names, and every job uses the standard preamble and `logs/`
  - every Markdown link resolves
  - the migration script works
  - every environment a job activates has a container
- **Container route checked by hand** with all the real images and a small synthetic 3-genome
  project. These ran successfully:
  - step 01: synteny, species check, contig redundancy, contig size profile
  - step 02: RepeatMasker, TE-Aid, the TE filtering tools, MCHelper Step 3
  - step 04: protein extraction
  - step 06: link + get_pangenes

## 8. Not yet verified / open items

- **Nothing on this branch has been run under real SLURM** or with the cluster's own apptainer.
  The tests and the container checks ran on a workstation.
- Containers: anything outside the project, the pipeline and `$HOME` (databases, the BUSCO
  lineage, genome symlink targets) must be listed in `CONTAINER_BINDS`.
- Step 06: the `<input>_pangenes/` output location was inferred from the Tgondii report and
  matches the local get_pangenes run, but hasn't been checked on the cluster.
- Tgondii's dropped-genes summary script was never committed. The job skips it with a warning.
- Step 03, InterProScan and BUSCO still depend on the cluster's modules and shared containers.
- `results/` was left as it is, and is planned to be removed.
- Building `mchelper.sif` needs `--fakeroot` (or root). If the cluster allows neither, build it
  on a Linux machine and copy it over.
