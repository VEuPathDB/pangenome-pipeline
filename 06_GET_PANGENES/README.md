# 06 – GET_PANGENES pipeline

***Purpose***: Generate the pangenes from the annotated genomes, i.e. the final GFFs of the
[basket annotation (step 05)](../05_Basket_annotation_integration/README.md).

From the [GET_PANGENES paper](https://link.springer.com/article/10.1186/s13059-023-03071-z):

"The protocol get_pangenes performs whole genome alignments (WGA) to call syntenic gene models based on coordinate overlaps. A benchmark with small and large plant genomes shows that pangenes recapitulate phylogeny-based orthologies and produce complete soft-core gene sets."

Source code and documentation: <https://github.com/Ensembl/plant-scripts/tree/master/pangenes>.

The associated publication is:

Contreras-Moreira et al. 2023. **GET_PANGENES: calling pangenes from plant genome
alignments confirms presence-absence variation.** *Genome Biology* 24, 223.
<https://doi.org/10.1186/s13059-023-03071-z>

## Table of contents

- [Installing get_pangenes](#installing-get_pangenes)
- [The two modes](#the-two-modes)
- [Step 1: put each genome's GFF and FASTA side by side](#step-1-put-each-genomes-gff-and-fasta-side-by-side)
- [Step 2: run get_pangenes](#step-2-run-get_pangenes)
- [Both steps in one go](#both-steps-in-one-go)
- [Checking it worked](#checking-it-worked)

## Installing get_pangenes

[GET_PANGENES](https://github.com/Ensembl/plant-scripts/tree/master/pangenes) can be installed via `micromamba`
(an environment called `get_pangenes` here). Step 2 reads the environment's `bin/` folder from `GET_PANGENES_ENV`
at the top of its script: edit that line to point to your own.

## The two modes

Step 05 writes two versions of every genome's annotation: `keep_similar` (a genome keeps its own model when it is
similar to the panel model) and `prefer_panel_model` (the panel model whenever it works). This step runs **once per
mode**. Both scripts choose the mode with the variable `MODE`, given when the job is submitted:

```bash
sbatch --export=ALL,MODE=prefer_panel_model <script>
```

`--export=ALL` passes your terminal's environment to the job (as sbatch does by default) and `MODE=...` adds the
mode. Without `MODE` the scripts use `keep_similar`. Each mode has its own folders, so the two never touch each other.

## Step 1: put each genome's GFF and FASTA side by side

Script: [`scripts/link_for_get_pangenes_basket.slurm.sh`](scripts/link_for_get_pangenes_basket.slurm.sh)
(1 CPU, 1 GB, a few minutes).

`get_pangenes` pairs a genome's GFF and FASTA by name: both must sit in the same folder with the same name, only the
extension differs (`.gff` / `.fna`). The step 05 GFFs are otherwise ready (feature type `gene`, IDs unique across the
panel). The script copies each GFF into the folder and links the genome FASTA next to it. A gene on a contig that is
not in that FASTA (e.g. a contig removed as redundant) is left out, and the log counts them; it should be 0 when step
05 and this step use the same genomes.

Settings, at the top of the script:

| Variable | Default | Meaning |
|---|---|---|
| `MODE` | `keep_similar` | the step 05 mode (set it with `--export`, see above) |
| `ANNOTATION_DIR` | `basket_annotation/gff_${MODE}_final` | the final GFFs of step 05 |
| `ANNOTATION_SUFFIX` | `_basket_annotation.gff` | GFF name = `<genome><suffix>` |
| `GENOME_DIR` | `genomes` | the genome FASTA files; must be the same genomes as in step 05 (*T. gondii*: `genomes_nonredundant`) |
| `GENOME_SUFFIX` | `.fna` | FASTA name = `<genome><suffix>` |
| `ORDER` | `genome_order.txt` | the genome list used throughout the pipeline |
| `OUT_DIR` | `for_get_pangenes_basket_${MODE}` | where the pairs are written (step 2's input) |

```bash
sbatch --export=ALL,MODE=prefer_panel_model 06_GET_PANGENES/scripts/link_for_get_pangenes_basket.slurm.sh
```

## Step 2: run get_pangenes

Script: [`scripts/run_get_pangenes_basket_job.sh`](scripts/run_get_pangenes_basket_job.sh) (16 CPUs, 64 GB, up to
48 h). It needs step 1 to have finished for the same mode.

The output goes to `<INPUT_DIR>_pangenes/` (e.g. `for_get_pangenes_basket_prefer_panel_model_pangenes/`). **If that
folder already exists, get_pangenes re-uses the files cached in it**: rename it before a new run, e.g.

```bash
mv for_get_pangenes_basket_prefer_panel_model_pangenes for_get_pangenes_basket_prefer_panel_model_pangenes_previous
```

Settings, at the top of the script:

| Variable | Default | Meaning |
|---|---|---|
| `MODE` | `keep_similar` | the step 05 mode (set it with `--export`) |
| `GET_PANGENES_ENV` | `/home/aludwig/hc-storage/micromamba/envs/get_pangenes/bin` | `bin/` of the get_pangenes environment |
| `INPUT_DIR` | `./for_get_pangenes_basket_${MODE}` | step 1's `OUT_DIR` |
| `N_CPUS` | `16` | must match `--cpus-per-task` in the `#SBATCH` lines |
| `TAXA_THRESHOLD` | `0` | get_pangenes' `-t`: minimum number of genomes a cluster must be in to be reported; `0` reports all, including genome-private clusters |
| `TAXON_REGEX` | unset | get_pangenes' `-s` (optional, exported before submitting): a regex grouping genome names into taxa |

The script puts the environment's `bin/` on the `PATH` (get_pangenes.pl, its bundled Python and bedtools) and clears
`PERL5LIB`-style variables first, so Perl libraries from other environments are not mixed in.

```bash
sbatch --export=ALL,MODE=prefer_panel_model 06_GET_PANGENES/scripts/run_get_pangenes_basket_job.sh
```

## Both steps in one go

Submitting both at once would start step 2 before step 1 has written its folder. To chain them, submit step 2 with a
dependency on step 1 (it starts only when step 1 has finished without error):

```bash
MODE=prefer_panel_model
jid=$(sbatch --parsable --export=ALL,MODE=$MODE 06_GET_PANGENES/scripts/link_for_get_pangenes_basket.slurm.sh)
sbatch --dependency=afterok:$jid --export=ALL,MODE=$MODE 06_GET_PANGENES/scripts/run_get_pangenes_basket_job.sh
```

## Checking it worked

- **Step 1** (`link_for_get_pangenes_basket_<jobid>.out`): one line per genome, `<genome>: N genes; left out (contig
  not in the FASTA): 0`, then `15 genome(s) prepared in for_get_pangenes_basket_<mode>`. A line `ERROR: missing ...`
  names the file not found.
- **Step 2** (`pangenes_basket_<jobid>.log` / `.err`): ends with `get_pangenes completed (mode <mode>). Outputs
  (pangene_matrix*.tab, etc.) are in for_get_pangenes_basket_<mode>_pangenes`. `ERROR: ... not found - run
  link_for_get_pangenes_basket.slurm.sh first` means step 1 had not run (or not finished) for that mode.
- The pangene matrices (`pangene_matrix.tr.tab`, `pangene_matrix_genes.tr.tab`) are in a subfolder of the output
  folder; [step 07](../07_Post_GET_PANGENES_adjusting/README.md)'s jobs find them by name.

[Back to top level](../README.md)
