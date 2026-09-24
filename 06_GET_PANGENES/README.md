# 06 – GET_PANGENES pipeline

***Purpose***: Generate the pangenes from the annotated genomes.

From the [GET_PANGENES paper](https://link.springer.com/article/10.1186/s13059-023-03071-z):

"The protocol get_pangenes performs whole genome alignments (WGA) to call syntenic gene models based on coordinate overlaps. A benchmark with small and large plant genomes shows that pangenes recapitulate phylogeny-based orthologies and produce complete soft-core gene sets."

Source code and documentation: <https://github.com/Ensembl/plant-scripts/tree/master/pangenes>.

The associated publication is:

Contreras-Moreira et al. 2023. **GET_PANGENES: calling pangenes from plant genome
alignments confirms presence-absence variation.** *Genome Biology* 24, 223.
<https://doi.org/10.1186/s13059-023-03071-z>

## Installing get_pangenes

[GET_PANGENES](https://github.com/Ensembl/plant-scripts/tree/master/pangenes) can be
installed via `micromamba`.

## Table of contents

1. [Link each genome's GFF and FASTA into one directory](#step-1-link-each-genomes-gff-and-fasta-into-one-directory)
2. [Run get_pangenes](#step-2-run-get_pangenes)

Two stages, in order: get every genome's GFF and FASTA into one directory, paired up by
filename the way `get_pangenes` expects, then run `get_pangenes` itself.

> The relabeling this step used to do (`protein_coding_gene`/`pseudogene` → feature type
> `gene`, `pseudogenic_transcript` → `transcript`) is no longer needed here - step 05's own
> [7a](../05_Cross-genome_annotation_integration/README.md#7a-filter-drop-unsupported-genes-before-get_pangenes)/[7b](../05_Cross-genome_annotation_integration/README.md#7b-no-filtering-tag-decide-after-get_pangenes)
> already write their output in exactly that shape (and already renumber every gene ID to be
> unique across the whole panel, which this step never did). Non-coding RNA genes are also
> already left correctly untouched by 7a/7b - typed away from `gene` back in
> [step 05's step 3](../05_Cross-genome_annotation_integration/README.md#step-3-prepare-the-gffs-for-liftoff),
> so `get_pangenes` never sees them as candidates, without needing any `source`-column retag.


## Dependencies

| Program | Used by | Suggestion |
|---|---|---|
| bash (`ln`) | Step 1 | nothing to install |
| get_pangenes (with its bundled Perl/Python and bedtools) | Step 2 | micromamba env (`get_pangenes` here), see [Installing get_pangenes](#installing-get_pangenes) |

The Step 2 scripts read the environment's `bin/` directory from `GET_PANGENES_ENV` at the top of
each script - edit that line to point to your own.

## Steps

### Step 1: Link each genome's GFF and FASTA into one directory

Input: each genome's `<genome>_final_annotation_filtered.gff` (from
[step 05's 7a](../05_Cross-genome_annotation_integration/README.md#7a-filter-drop-unsupported-genes-before-get_pangenes))
or `<genome>_final_annotation_tagged.gff` (from
[7b](../05_Cross-genome_annotation_integration/README.md#7b-no-filtering-tag-decide-after-get_pangenes)),
whichever approach you're running with.

`get_pangenes` pairs up a genome's GFF and FASTA by filename - both need to sit in the same
directory with the **exact same basename**, differing only in extension (`.gff` vs `.fna`).
7a/7b's output is already `get_pangenes`-ready in every other respect (feature types, unique
IDs), so all that's actually left is this pairing - nothing to relabel or rewrite any more.
This step just symlinks both into place (no file is copied), one job, instant.

**Separate scripts for 7a and 7b** - both can run independently, each into its own directory,
so nothing stops you running `get_pangenes` on both and comparing the two pangene calls:
[`scripts/link_for_get_pangenes_7a.slurm.sh`](scripts/link_for_get_pangenes_7a.slurm.sh) ·
[`scripts/link_for_get_pangenes_7b.slurm.sh`](scripts/link_for_get_pangenes_7b.slurm.sh).

#### Dependencies

None - plain bash, no AGAT/Python/etc.

#### What you need before starting

- Every genome's own `<genome>_final_annotation_filtered.gff` (7a) or `_tagged.gff` (7b), from
  step 05's 7a/7b `OUT_DIR`.
- Every genome's own genome FASTA, in `genomes/` (step 01's assemblies folder).
- **The same genome order** used throughout this pipeline (`genome_order.txt`, first = reference).

#### Running it

**Edit the settings at the top of whichever script matches what you're running** (both have the
same variables, pre-filled with each approach's own defaults):

| Variable | 7a default | 7b default |
|---|---|---|
| `ANNOTATION_DIR` | `filtered_annotation` | `tagged_annotation` |
| `ANNOTATION_SUFFIX` | `_final_annotation_filtered.gff` | `_final_annotation_tagged.gff` |
| `GENOME_DIR` | `genomes` | `genomes` |
| `GENOME_SUFFIX` | `.fna` | `.fna` |
| `ORDER` | `genome_order.txt` | `genome_order.txt` |
| `OUT_DIR` | `for_get_pangenes_7a` | `for_get_pangenes_7b` |

**Run** (either or both - they don't touch each other's files):
```bash
sbatch 06_GET_PANGENES/scripts/link_for_get_pangenes_7a.slurm.sh
sbatch 06_GET_PANGENES/scripts/link_for_get_pangenes_7b.slurm.sh
```

#### Checking it worked

- `<OUT_DIR>/<genome>_for_get_pangenes.gff` and `.fna` (both symlinks) should exist for every
  genome once the job finishes.
- The job's log ends with `N genome(s) linked into <OUT_DIR>. Outputs are in <OUT_DIR>` or an
  `ERROR: missing ...` line naming exactly which file wasn't found.

### Step 2: Run get_pangenes

Input: each genome's `<genome>_for_get_pangenes.gff` + matching `.fna`, from step 1 above.

Same split as step 1 - each variant reads its own `INPUT_DIR` and never touches the other's:
[`scripts/run_get_pangenes_7a_job.sh`](scripts/run_get_pangenes_7a_job.sh) ·
[`scripts/run_get_pangenes_7b_job.sh`](scripts/run_get_pangenes_7b_job.sh).

#### Dependencies

- **get_pangenes**, from the micromamba environment set in `GET_PANGENES_ENV`. The script puts
  that environment's `bin/` on the `PATH` and uses its `get_pangenes.pl`, its bundled Python and
  its `bedtools`. It also clears any `PERL5LIB`-style variables first so Perl libraries from other
  environments do not get mixed in.

#### What you need before starting

- Every genome's `<genome>_for_get_pangenes.gff` + `.fna` in `for_get_pangenes_7a`/`_7b` -
  step 1's output.
- `get_pangenes` itself installed (see "Installing get_pangenes" above).

#### Running it

**Edit the settings at the top of whichever script matches what you're running** (same
variables in both, pre-filled with each approach's own `INPUT_DIR`):

| Variable | Meaning |
|---|---|
| `GET_PANGENES_ENV` | `bin/` directory of the micromamba environment from the install step above |
| `INPUT_DIR` | Directory with every genome's `<genome>_for_get_pangenes.gff` (+ matching `.fna`) - `./for_get_pangenes_7a`/`_7b` by default |
| `N_CPUS` | Should match `--cpus-per-task` in the `#SBATCH` header above it |
| `TAXA_THRESHOLD` | `get_pangenes`' `-t` - minimum number of genomes a pangene cluster must appear in to be reported at all; `0` reports everything, including genome-private (singleton) clusters |
| `TAXON_REGEX` (optional, exported before submitting) | `get_pangenes`' `-s` - a regex to group genome names into taxa; if unset, every genome is treated as its own taxon |

**Run** (either or both):
```bash
sbatch 06_GET_PANGENES/scripts/run_get_pangenes_7a_job.sh
sbatch 06_GET_PANGENES/scripts/run_get_pangenes_7b_job.sh
```

#### Checking it worked

- The job's log (`pangenes_7a_<jobid>.log`/`.err` or `pangenes_7b_...`) should end with
  `get_pangenes (7a/7b) completed. Outputs (pangene_matrix*.tab, etc.) are in <INPUT_DIR>`.
- `get_pangenes` writes its output (`pangene_matrix*.tab` and its other result files) directly
  into `INPUT_DIR`, alongside the GFFs/FASTAs it read - not a separate output directory.

[Back to top level](../README.md)
