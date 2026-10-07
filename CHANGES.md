# Changes on the `add-python-tooling` branch

What's different from `main`, for anyone testing the branch. The folder layout is **the same as
`main`**: step folders `NN_<Step>/scripts/`, the README in each step folder, the pipeline cloned
into the species project, and every job run from the project folder and writing where it did on
`main`. Each commit message has the full detail (`git log main..add-python-tooling`).

## Quick start for testing

```bash
git clone git@github.com:adriludwig/pangenome-pipeline.git <project>   # or git pull in an existing clone
cd <project>
git switch add-python-tooling

bash setup_project.sh                              # genomes/, gff_originals/ and the jobs' log folders
cp pipeline_config.example.sh pipeline_config.sh   # fill in what you use
uv sync --frozen --no-dev                          # the pipeline's Python (install uv first: see README)
sbatch --export=ALL,ENV_BACKEND=path scripts/test_python_env.slurm.sh   # log ends with "ALL OK"
```

## 1. Configuration instead of hardcoded values

**`pipeline_config.sh`** (at the project root, from `pipeline_config.example.sh`, not committed)
holds every cluster-specific path that used to be written into the scripts. Any setting can also
be given per job with `sbatch --export=ALL,NAME=value`, and a job that needs an empty setting
stops with a message naming it.

**Species defaults removed.** These were silently set to one species' values. They now have to be
set, and the job stops if they aren't:

| Setting | Used by | Was |
|---|---|---|
| `SPECIES`, `STRAIN` | step 03 annotate / predict / train | A. fumigatus Af293 |
| `SOURCE_GENOME`, `ORGANELLE_CONTIGS` | step 01 organelle screen | AF293's mitochondrion |
| `LABEL` | steps 07 and 08 jobs | `Toxo7b` |
| `ID_PREFIX` | step 05 `tag_final_annotation_job.sh` | `Afu` |

The `#SBATCH` partition/account can be overridden for every job with
`export SBATCH_PARTITION=... SBATCH_ACCOUNT=...` before `sbatch`.

## 2. Tool environments: micromamba is optional

Jobs call `activate_env <name>`, which uses, per environment, the first that is set up: a
container in `$CONTAINER_DIR`, micromamba, conda, or whatever is already on `PATH`. The job's log
says which (`activate_env synteny: container`).

- **Python** runs from a uv venv pinned by `uv.lock` (`uv sync --frozen --no-dev`), so every
  cluster gets the same package versions. Tested on a compute node with
  `scripts/test_python_env.slurm.sh`.
- **`containers/pull_containers.sh DIR`** downloads pinned images (biocontainers, `dfam/tetools`)
  for every environment; **`containers/mchelper.def`** builds MCHelper's, which has no published
  image. Folders such as `synteny/` hold one image per tool, used together.
- The old settings still work: `GET_PANGENES_ENV`, `TE_AID`, `MCHELPER_PY` and micromamba are
  used when set.
- Not covered: step 03 (funannotate), InterProScan and BUSCO still use the cluster's shared
  containers and modules.

## 3. Python code

- **Shared helpers in `lib/pangenome_common/`** instead of a copy per script: GFF finding and
  attribute parsing, the `ID_RE` / `PARENT_RE` / … patterns, biotypes, intervals, genome and ID
  lists, PAF reading, `load_extra_genes`. Step 05 shares `annotation_common.py`, Get-PanLoci
  `panloci_common.py`. Scripts find `lib/` from their own location.
- Byte-identical duplicate scripts removed (Get-PanLoci); the four step 06 wrappers share
  `get_pangenes_common.sh`; the Afum-only RepeatMasker fix and Tgondii scripts moved to
  `species_specific/`.
- Formatted and linted (ruff, isort, mypy).

## 4. Changes that can change results

Bug fixes, but outputs can differ from a run on `main`:

| Step | Change |
|---|---|
| 01 | `contig_redundancy.slurm.sh` never passed `--ref`, so for any reference other than AF293 every contig showed 0 % placement with no error. It now uses the first genome in `genome_order.txt`. |
| 02 | Gene descriptions decode `+` as a space (`hypothetical+protein` → `hypothetical protein`), as steps 01 and 05 already did. |
| 01, 02 | The description is the first *non-empty* description attribute; trailing whitespace / `\r` is stripped from attribute values. |
| 03 | Jobs submitted from `03_denovo_Genome_Annotation/` (as its README says) now find `lib/` above it; before, every funannotate job stopped at that line. |
| 05 | Gene renumbering broke ties (same start position, common for Liftoff extra copies) in hash order, so IDs could swap between runs; ties now go by the original gene ID. |
| 05 | `ID=`, `Parent=`, `extra_copy_number=`, `valid_ORF=` match whole attributes only (`sequence_ID=` was read as the ID if it came first). |
| 05 | Liftoff `coverage=1e-05` was read as `1.0` and passed `--min-rescued-coverage`; it's now parsed in full. |

## 5. Tooling and tests

- `pyproject.toml` + `uv.lock` pin the dependencies (a test fails if the lock goes stale);
  `nox -s tests` runs ruff, format, isort, mypy, toml checks and pytest.
- **284 tests**: every script compiles and parses its CLI, shell syntax, unit tests for the shared
  helpers, step 05 end-to-end snapshots, jobs run from a fake project, and guards against
  hardcoded paths, references to missing scripts, and environments without a container.

## 6. Not yet verified / open items

- Only the Python test job has run under real SLURM. The container route was checked by hand on a
  workstation (steps 01, 02, 04, 06 with the real images), not on the cluster.
- Containers: anything outside the project folder and the pipeline (databases, the BUSCO lineage,
  genome symlink targets) must be listed in `CONTAINER_BINDS`.
- The `low` partition can pre-empt and requeue jobs; long jobs restart from scratch.
- Tgondii's dropped-genes summary script was never committed (as on `main`).
- Building `mchelper.sif` needs `--fakeroot` (or root). If the cluster allows neither, build it on
  a Linux machine and copy it over.
