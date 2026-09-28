# 08 – Pangenome Visualisation

***Purpose***: Turn the `get_pangenes` result into a fixed set of figures, so every species and every run
is looked at the same way. The figures answer two questions: what the pangenome looks like, and how much of
its accessory part (the cloud) may reflect real biology and how much may come from assembly or annotation problems.

## Table of contents

1. [Dependencies](#dependencies)
2. [Input](#input)
3. [InterPro entry list](#interpro-entry-list)
4. [Plots and their purpose](#plots-and-their-purpose)
5. [Running it](#running-it)
6. [Adapting the scripts](#adapting-the-scripts)
7. [Outputs](#outputs)
8. [Additional options](#additional-options)

## Dependencies

| Program | Used by |
|---|---|
| python3 (3.9 or newer) + numpy + matplotlib | every `plot_*.py` script |
| scipy (1.11 or newer) | `plot_pancore.py` (curve fits), `plot_domain_enrichment.py` (Fisher test, p-value correction), `plot_occupancy_properties.py` and `plot_occupancy_contigs.py` (Tukey test) |
| pandas + seaborn | `plot_pocs.py` (POCS heatmap) |

The jobs are hardcoded with the author's own micromamba environment (paths and env name `synteny`) - edit the
`### ADAPT ###` block at the top of each job to point to your own, see [Adapting the scripts](#adapting-the-scripts).
Easiest fix: create an env with the same name and every package above installed:

    micromamba create -n synteny -c conda-forge python numpy scipy matplotlib pandas seaborn

If `synteny` already exists from [step 01](../01_Genome_Selection_QC/README.md#dependencies), it has most of them; add what is
missing with `micromamba install -n synteny -c conda-forge scipy seaborn`.


### InterPro entry list

The enrichment plot names each InterPro entry by its IPR accession. The names come from the InterPro file `entry.list`,
a tab-separated table with the accession, the type and the name of every entry. Only `plot_domain_enrichment.py` reads it.

Download it once and save it in a folder that all your
projects can read, and use the path of that folder in the command:

```bash
curl -L -o /path/to/interpro_entry.list https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list
```

The file comes from the current InterPro release, which changes over time, so note the download date.

**Tell the scripts where it is.** Use one of these:

- Set `INTERPRO_ENTRY_LIST` in the project's `pipeline_config.sh` (see the
  [main README](../README.md#getting-the-pipeline-into-a-project)) - both jobs use it by default:

  ```bash
  INTERPRO_ENTRY_LIST="/path/to/interpro_entry.list"
  ```

- On the command line of `plot_all_job.sh`: `--entry-list /path/to/interpro_entry.list`, or the
  `ENTRY_LIST=` line of [`plot_domain_enrichment_job.sh`](scripts/plot_domain_enrichment_job.sh).

With neither, the jobs warn and label terms by IPR id.

- Running `plot_domain_enrichment.py` by hand: `--entry-list /path/to/interpro_entry.list`. Without it, the script
  downloads the file next to the figure, as `<out-prefix>_entry.list`.

If the file is not at that path, the two jobs try to download it there. If that fails as well, the plot still runs and
labels the entries by their IPR accession.

## Input

| Input | From | Used by |
| :--- | :--- | :--- |
| `pangene_matrix.tr.tab` (rows = pangenes, columns = genomes, cells = number of genes) | `get_pangenes` | pan/core, occupancy, singletons |
| `pangene_matrix_genes.tr.tab` (same, cells = gene IDs) | `get_pangenes` | every plot that looks at genes |
| `POCS.matrix.tab` | `get_pangenes`, same folder as the matrices | POCS heatmap |
| the GFFs given to `get_pangenes` (`<genome>_for_get_pangenes.gff`) | [step 06](../06_GET_PANGENES/README.md), gene IDs must match the matrix | domains, protein length, `valid_ORFs`, gene source, `intersect=` |
| `contig_sizes.tsv` (genome, contig, length) | [step 01](../01_Genome_Selection_QC/README.md), contig size profile | contig size plot |
| InterPro `entry.list` (IPR id → name) | downloaded once from EBI | names in the enrichment plot |

`get_pangenes` writes the matrices in a subfolder whose name depends on the genome and the options it was run
with, so the job looks for them by file name. There are two sets of pangenes, and the same plots are made for each:

- **7a**, with filtering: the GFFs were filtered before `get_pangenes` ([step 05](../05_Cross-genome_annotation_integration/README.md#7a-filter-drop-unsupported-genes-before-get_pangenes)).
- **7b**, no filtering ([step 05](../05_Cross-genome_annotation_integration/README.md#7b-no-filtering-tag-decide-after-get_pangenes)): the GFFs were only tagged, nothing removed. It carries the `intersect=` tag, which one plot needs.

**Occupancy categories.** Occupancy is the number of genomes (N) a pangene is present in. The categories are cut in
numbers of genomes, not in percentages. Each pangene falls in exactly one, decided in this order (`categorize` in
[`plot_occupancy.py`](scripts/plot_occupancy.py), used by every plot that has categories):

| Category | Default rule | Parameter |
| :--- | :--- | :--- |
| Core | present in all N genomes | `--core-missing 0`: missing from at most 0 genomes |
| Soft-core | not core, present in N-1 genomes | `--softcore-missing 1`: missing from at most 1 genome |
| Cloud | not core or soft-core, present in 1 or 2 genomes | `--cloud-max 2`: present in at most 2 genomes |
| Shell | everything else | |

What the defaults give for some panel sizes:

| Genomes (N) | Cloud | Shell | Soft-core | Core |
| ---: | :---: | :---: | :---: | :---: |
| 15 | 1-2 | 3-13 | 14 | 15 |
| 17 | 1-2 | 3-15 | 16 | 17 |
| 30 | 1-2 | 3-28 | 29 | 30 |

The defaults do not depend on N. With many genomes they may not suit the panel, and can be adapted with the three
parameters (see [Running it](#running-it)). The plots are always drawn with cloud on the left and core on the right.

## Plots and their purpose

Files start with the label given to the run (for example `Afum7b`). Every figure also comes as `.svg`.

| Plot | File | What it shows | Purpose |
| :--- | :--- | :--- | :--- |
| Pangenome and core genome size | `_pancore`, `_pancore_separate` | Genomes are added one at a time in 100 random orders. The pangenome and the core size are drawn at each step, with the fitted curves (Tettelin for the pan, Tettelin and Willenbrock for the core). The `.log` has the fit parameters. | Is the pangenome open (still growing with each genome) and is the core stable, or is the panel too small to tell? |
| Pangene occupancy | `_occupancy` | One bar per occupancy (1 to N genomes), coloured by category. The categories are decided as in the table above. | The overall composition: how much is core, soft-core, shell and cloud. |
| Singletons per genome | `_singletons` | Pangenes found in only one genome, counted per genome. | Genome-specific content, and outliers. A genome with many more singletons than the others may be divergent, or carry redundant contigs, contamination or poor annotation. When two or more genomes are very close, what they share is present in more than one genome, so their singletons are automatically fewer. |
| POCS heatmap | `_POCS_heatmap` | Percentage of conserved sequences between every pair of genomes, clustered. | How related the genomes are, and which one sits apart. |
| Gene properties by occupancy | `_occupancy_properties`, `_occupancy_properties_log` | Top: share of genes with an InterPro/PFAM domain. Bottom: protein length. One value per pangene, by occupancy. | Short genes with no domain tend to include annotation artefacts and gene fragments, but they can also be real genes. The plot shows how much the cloud differs from the core in these properties. |
| Contig size by occupancy | `_contig_classes` | For each category, the share of its pangenes in each contig size class (median contig of the pangene's genes). | A cloud that sits on short contigs may reflect leftovers of the assembly (redundant or contaminant contigs), although short contigs can also carry real accessory genes. |
| InterPro term enrichment | `_domain_enrichment` | Which InterPro entries are found more often in the core or in the cloud than in the other pangenes (one-sided Fisher test, Benjamini-Hochberg adjusted). Dot size is the gene ratio, colour the adjusted p-value. Entries carried by fewer than 5 pangenes in the whole universe are not tested, and an entry is drawn only if at least 3 pangenes of the category carry it. | What functions appear to make up the core and the cloud. An entry that dominates the cloud can point to an artefact, for example a very long gene annotated in pieces or an organelle gene, and is worth checking. |
| Invalid ORFs by occupancy | `_occupancy_invalid_only` (count), `_occupancy_invalid_percent` (percentage), `_invalid_share` | Pangenes with more than 50 % of invalid genes, by occupancy. `_invalid_share` has one bar per category, from cloud to core, split by the share of invalid genes in each pangene: none, up to half, more than half but not all, all. A gene is invalid when Liftoff found no valid ORF (`valid_ORFs=0`) or the source annotation calls it a pseudogene. | Shows how much of each category is made of genes with a broken ORF or a pseudogene mark. Invalid pangenes concentrated in the cloud suggest broken transfers or pseudogenes rather than functional accessory genes, although a pseudogene can also be a real, recent loss. |
| Unsupported overlapping genes by occupancy | `_occupancy_unsupported_only` (count), `_occupancy_unsupported_percent` (percentage), `_unsupported_share` | Pangenes with more than 50 % of genes that are not from the reference layer (`main`), overlap another gene and have neither a domain nor a specific description. **7b only.** | Shows possible false genes that may exist only because they overlap a real one. If many are in the core, they may be conserved only through the gene they overlap. |

Notes on how the pangenes are counted:

- **Gene properties** and **contig size** use only pangenes with exactly one gene in each genome they are present in.
  A pangene with several genes in a genome has no single protein length or contig, and it may hold distinct genes.
  The job log says how many pangenes were left out.
- **Invalid ORFs** and **unsupported overlapping genes** count every gene of every pangene, several genes in one
  genome included. Genes taken from a genome's own annotation have no ORF result and can only be marked by the
  pseudogene biotype; pangenes with no gene that can be judged are left out.
- The Tukey p-value written on the gene-properties plot treats the genes of one pangene as independent samples,
  although they come from different genomes, so it is optimistic.
- The enrichment counts only pangenes with at least one InterPro entry. A pangene carries an entry if at least one of its
  genes, in any genome, has it; how many genes have it does not matter. A core pangene collects entries from the genes of
  every genome, which may favour the core side.

## Running it

On the cluster, one job makes all the plots. Run it once for each set of pangenes (7a and 7b), from the project folder.
The log folder must exist before submitting:

```bash
mkdir -p pangenome_plots/logs
sbatch 08_Pangenome_Visualisation/scripts/plot_all_job.sh --set 7b --label Afum7b
sbatch 08_Pangenome_Visualisation/scripts/plot_all_job.sh --set 7a --label Afum7a
```

The two commands differ in two things:

- **The folders they read.** `--set 7b` reads `for_get_pangenes_7b_pangenes` and `for_get_pangenes_7b`. `--set 7a` reads
  `for_get_pangenes_7a_pangenes` and `for_get_pangenes_7a`.
- **The plots they make.** `--set 7b` makes every plot, including the unsupported overlapping genes plots, which need the
  `intersect=` tag that only the 7b GFFs have. `--set 7a` makes all the others.

`--label` only names the output files (`Afum7b_...` and `Afum7a_...`), so the two runs do not overwrite each other.

**Folders with other names.** If your folders are not called `for_get_pangenes_<set>_pangenes` and
`for_get_pangenes_<set>`, do not use `--set`. Give the folders yourself, and for the 7b set also add `--overlap yes`:

```bash
sbatch 08_Pangenome_Visualisation/scripts/plot_all_job.sh --pangenes-dir <get_pangenes folder> --gff-dir <GFF folder> --label Afum7b --overlap yes
```

For a 7a set do the same without `--overlap`. `--overlap` is only needed in this case; with `--set` it is decided for you.

[`plot_all_job.sh`](scripts/plot_all_job.sh) options (`--help` lists them all):

| Option | Meaning |
| :--- | :--- |
| `--label` | Prefix of every output file. Required. |
| `--set 7a` / `--set 7b` | Shortcut that chooses the folders (`--pangenes-dir for_get_pangenes_<set>_pangenes --gff-dir for_get_pangenes_<set>`). `7b` also makes the unsupported overlapping genes plots. |
| `--pangenes-dir`, `--gff-dir` | The `get_pangenes` result folder and the GFFs given to it, when the names are different. |
| `--matrix-subdir` | Only when the result folder holds more than one `pangene_matrix.tr.tab`: the subfolder to use. |
| `--gff-suffix` | GFF name is `<genome><suffix>` (default `_for_get_pangenes.gff`). |
| `--overlap yes` | Makes the unsupported overlapping genes plots. Only for the 7b set, when you give `--pangenes-dir` and `--gff-dir` yourself instead of `--set 7b`. Default: not made. |
| `--contig-sizes` | Path of `contig_sizes.tsv` (default `01_Genome_Selection_QC/contig_size_profile/contig_sizes.tsv`). |
| `--entry-list` | Path of the InterPro `entry.list`; downloaded there if missing. Download it once on the login node, because compute nodes may have no internet. |
| `--core-missing`, `--softcore-missing`, `--cloud-max` | How the occupancy categories are cut, in numbers of genomes (default 0, 1, 2; see above). Applied to every plot with categories. |
| `--categories`, `--alpha`, `--min-total`, `--min-count`, `--top` | Enrichment settings (default `Core,Cloud`, 0.05, 5, 3, 15): the categories tested, the adjusted p-value cut, the minimum size of an entry in the whole universe to be tested, the minimum pangenes of the category carrying an entry to draw it, and the entries drawn per category. |
| `--out-dir` | Output folder (default `pangenome_plots`). |

The job looks for `pangene_matrix.tr.tab` by name below `--pangenes-dir`. If it finds none, or more than one, it
stops and says so. If one plot fails, the others still run; the job lists the failed ones at the end and exits with an
error. Set `SKIP_ENV=1` to skip the micromamba block when already inside an environment.

To make only the first four plots (pan/core, occupancy, singletons, POCS) on a laptop, in seconds (the category cuts can be
set with the environment variables `CORE_MISSING`, `SOFTCORE_MISSING` and `CLOUD_MAX`):

```bash
bash 08_Pangenome_Visualisation/scripts/plot_all.sh -i <folder>/pangene_matrix.tr.tab -o <out_dir> -l <label>
```

## Adapting the scripts

The Python scripts need no change. What is specific to one cluster and one user is in the jobs.

In [`plot_all_job.sh`](scripts/plot_all_job.sh):

- The `#SBATCH` lines at the top: partition, memory, time and the log folder (`pangenome_plots/logs`).
- The block between `### ADAPT ###` and `### END ADAPT ###` activates the environment named by `ENV_SYNTENY`; the
  micromamba paths come from the project's `pipeline_config.sh`. Set `SKIP_ENV=1` to skip it when already inside an
  environment.
- The InterPro `entry.list`: `INTERPRO_ENTRY_LIST` in `pipeline_config.sh` (see
  [InterPro entry list](#interpro-entry-list)), or pass your own path with `--entry-list`.
- The default of `--contig-sizes` (`01_Genome_Selection_QC/contig_size_profile/contig_sizes.tsv`), if the step 01 table
  is elsewhere.

The jobs of the single plots (`plot_domain_enrichment_job.sh`, `plot_occupancy_properties_job.sh`) have the same
`#SBATCH` lines and environment block, and their folders, label and settings are variables at the top of each file.

## Outputs

One `.png` and one `.svg` per figure, plus a `.tsv` with the values behind it (`.log` for the pan/core fits):

```
pangenome_plots/
├── <label>_pancore.{png,svg,log}   <label>_pancore_separate.{png,svg,log}
├── <label>_occupancy.{png,svg,tsv}
├── <label>_singletons.{png,svg,tsv}
├── <label>_POCS_heatmap.{png,svg}
├── <label>_occupancy_properties.{png,svg,tsv}   <label>_occupancy_properties_log.{png,svg,tsv}
├── <label>_contig_classes.{png,svg,tsv}
├── <label>_domain_enrichment.{png,svg,tsv}
├── <label>_occupancy_invalid_only.{png,svg,tsv}   <label>_occupancy_invalid_percent.{png,svg,tsv}   <label>_invalid_share.{png,svg,tsv}
└── <label>_occupancy_unsupported_only.{png,svg,tsv}   <label>_occupancy_unsupported_percent.{png,svg,tsv}   <label>_unsupported_share.{png,svg,tsv}   (7b only)
```

The `_pangenes.tsv` files next to the invalid and unsupported plots list every pangene with its share of bad genes.

## Additional options

The job makes one fixed version of each figure. Two of the scripts can also draw other versions of the same figure, or
change how the pangenes are counted. The job does not make these; to get one, run the script by hand, with the environment
active. For example, the contig size plot as boxplots, for the 7b set:

```bash
python3 08_Pangenome_Visualisation/scripts/plot_occupancy_contigs.py \
  --genes-matrix <get_pangenes folder>/pangene_matrix_genes.tr.tab \
  --gff-dir for_get_pangenes_7b --gff-suffix _for_get_pangenes.gff \
  --contig-sizes 01_Genome_Selection_QC/contig_size_profile/contig_sizes.tsv \
  --style boxplot -o pangenome_plots/Afum7b_occupancy_contigs
```

The other arguments are the same as in the job. What each option does:

**Contig size plot** (`plot_occupancy_contigs.py`)
- `--style boxplot`: instead of the bars by size class, one boxplot per occupancy with the median contig length of each
  pangene, and a Tukey p-value for core against cloud. `--yscale log` or `--yscale linear` sets the axis.

**Invalid ORFs and unsupported overlapping genes plots** (`plot_occupancy_invalid.py`; `--mark overlap` for the second)
- `--single-copy`: leaves out the pangenes with more than one gene in a genome.
- `--min-share`: the share of bad genes above which a pangene counts as "mostly" (default 0.5, half).

[↑ Back to top](#08--pangenome-visualisation) &nbsp;|&nbsp; [Back to top level](../README.md)
