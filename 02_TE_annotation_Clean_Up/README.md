# 02 – TE Annotation and Clean up

***Purpose:*** To identify, curate and annotate bona fide TEs across all genomes, preventing TE-derived sequences from biasing the pangenome analysis while retaining their annotation for further investigation. 

Dependencies are documented within each step below, in its own "Dependencies" note.

## Table of Contents

- [TE-related challenges in pangenome construction](#te-related-challenges-in-pangenome-construction)
- [TE library curation](#te-library-curation)
  - [1 - Repbase TE database preparation](#1---repbase-te-database-preparation)
  - [2 - TE Candidate Discovery Pipeline](#2---te-candidate-discovery-pipeline)
  - [3 - Manual Validation of the candidate TEs](#3---manual-validation-of-the-candidate-tes)
  - [4 - Complementary TE Screen pipeline](#4---complementary-te-screen-pipeline)
- [TE annotation](#te-annotation)
- [TE composition analysis](#te-composition-analysis)
- [Filtering TEs from gene annotations](#filtering-tes-from-gene-annotations)
  - [Quantifying TE Overlap Across Gene Models](#quantifying-te-overlap-across-gene-models)
  - [Removing TE-Derived Genes with a Kill List](#removing-te-derived-genes-with-a-kill-list)
- [Generating a TE-only GFF3](#generating-a-te-only-gff3)
- [Identifying Starships in fungal genomes](#identifying-starships-in-fungal-genomes)

---

## **TE-related challenges in pangenome construction**

Transposable Elements (TEs) are important components of eukaryote genomes. Only a few organisms so far are proved to be devoid of any kind of classic TEs, including Plasmodium falciparum. 

TEs need to be carefully identified before or after building the pangenome because, due to their ability to move around, they often appear as cloud genes. TE-related domains can also dominate enrichment results and mask patterns in the other genes. For this reason, TE-derived genes will be removed from the gene-focused pangenome analysis.

One challenge is that TEs may be inserted within introns without disrupting the coding sequence. Therefore, we should not exclude every gene that overlaps a TE, but only genes with substantial TE coverage within their CDS. Some gene models may also incorrectly incorporate TE sequences into the CDS, so these cases should be treated as potential annotation errors and examined more carefully when possible

One caveat is that TEs can occasionally be domesticated by the host genome and acquire a new cellular function. These genes may still be identified as TEs based on their protein domains, even if they are no longer mobile. Distinguishing domesticated TE-derived genes would require a closer examination of individual copies. For example, we could check whether they have lost the sequences required for recognition and mobilisation by the TE machinery (TIRs and TSDs); however, this might be difficult to do on a large scale. We could instead investigate known examples from the literature and examine genes containing TE-related domains that are conserved across all strains as their conservation for a long evolutionary time may suggest that they have acquired an important host function.

However, TEs should not be ignored. Although they may not introduce new protein functions into the pangenome, since their proteins mainly support their own movement, they can still have a major impact on the genome. For example, TE insertions can affect the expression of nearby genes and may help explain phenotypic differences among strains. Thus, the best approach is to construct a curated TE library for each species and eliminate them from the pangene matrix while keeping their annotation.

A related but opposite case is that of Starships, giant fungal mobile elements that carry a cargo of functional host genes. Because that cargo is made of real genes, Starships will not be removed from the pangenome; instead they will be identified so their cargo can be recognised as mobile-element-driven variation rather than independent gene gain/loss.

[↑ Back to top](#table-of-contents)

---

## **TE library curation**

Curating TE content is not a trivial task. Nevertheless, some level of curation is essential because RepeatModeler identifies repetitive sequence families *de novo*, and not every repetitive sequence it detects is necessarily a TE. Its output may also contain redundant families, host genes and other repetitive genomic regions that could be mistakenly classified as TEs.

We initially considered using [MycoMobilome](https://github.com/TobyBaril/MycoMobilome), a database of TE consensus sequences generated from thousands of publicly available fungal genomes using a standardised and reproducible workflow. Although the database includes several filtering steps and provides a useful resource for fungal TE annotation, it is described as an uncurated starting point. In our assessment, the library still contained considerable redundancy and some apparent host-gene sequences and lots of unconfirmed/non-TE sequences. We also identified a few examples of TE misclassification. Filtering and validating the entire library would probably require more work than constructing species-specific libraries from scratch. Using [Dfam](https://dfam.org/) alone resulted in an underestimation of TE content in our genomes. We therefore developed the methodology to identify and curate the TEs. 

For some species, however, curated TE sequences are already available in [Repbase](https://www.girinst.org/repbase/update/index.html). Repbase is a reference database of representative repetitive sequences from eukaryotic species, including consensus sequences for many known TE families and subfamilies. When available for the species being analysed, these sequences were used as a pre-curated starting library.

Some TEs may still be missing from our final libraries. We chose a conservative approach and retained only *bona fide* TEs that could be supported by characteristic protein domains, structural features, or both. Consequently, divergent or non-canonical TEs may have been missed. Non-autonomous elements, including SINEs and MITEs, are also likely to be underrepresented because they do not encode the protein domains commonly used to confirm autonomous TEs. Thus, the final libraries favour confidence over completeness: the retained sequences have strong support as TEs, but they may not represent the complete TE diversity of each genome. We can keep improving the TE libraries over time if we identify new TE sequences. 

### 1 - Repbase TE database preparation

For some species, a curated TE library is already available in [Repbase](https://www.girinst.org/repbase/update/index.html). In general, Repbase sequences are well curated and reliable. Thus, first step is to check for the presence of TEs in Repbase for the species of interest. However, the library might not present the complete TE repertoire of the species as some TE families may have been lost during evolution or may be specific to the strain(s) being analysed. Thus, if available, this library will be a starting point to get the final TE library with the pipeline described below. The idea is that you only need to proceed with the TE curation for a consensus sequence that is not yet in this pre-curated library.

The elements should be inspected to decide whether to include them in the pre-curated library. LTR elements are frequently represented split into LTR and internal region. In this case, the LTRs should be joined to the internal region to be added to the library. Moreover, LTR solos should not be included in the pre-curated library. The exclusion allows the pipeline to recover the full-length LTR elements from the genome. If no complete element carrying the particular LTR is retrieved, the LTR consensus sequence is not included in the final library.

#### Output: 
`Pre-curated TE library`

### 2 - TE Candidate Discovery Pipeline

Pipeline available here: [TE Candidate Discovery](https://github.com/adriludwig/TE-candidate-pipeline)

We developed a pipeline [TE Candidate Discovery](https://github.com/adriludwig/TE-candidate-pipeline) to automate the identification of high-confidence candidate TE families across **multiple genome assemblies** of the same species while minimising the amount of manual curation required. The final output is a FASTA library of candidate TEs that should be **manually inspected** and **validated** before being incorporated into a curated TE library for downstream analysis.

Analysing **multiple genomes** also provides a more comprehensive representation of the mobilome present in the species pangenome, increasing the probability of recovering TE families that may be absent, fragmented, or poorly assembled in any individual genome while avoiding redundant manual inspection of the same family recovered independently from different assemblies.

This approach is well suited for genome annotation. However, analyses that rely on precise sequence divergence, such as RepeatMasker landscape analyses used to infer the relative ages of TE insertions, may be influenced by this approach.

Moreover, the pipeline adopts a conservative strategy. Candidate sequences lacking detectable nucleotide or protein similarity to previously described TEs are excluded from the final library. This is because it is unlikely that a *bona fide* canonical TE would lack detectable similarity to any known TE at both the nucleotide and protein levels. Although some of these candidates may represent genuine non-canonical or non-autonomous TEs, demonstrating their transposable nature generally requires extensive manual curation.

#### Output: 

- `final_potential_new_TEs.fa` — the RepBase/TE-peptide-filtered candidate library, the one
  meant for validation and incorporation into the curated TE library.
- `all_MCHelper_curated_sequences_NR.cdhit80.fa` — the unfiltered non-redundant set, before
  that similarity filter. Worth checking separately if you're specifically interested in MITEs
  or solo LTRs: these non-autonomous/incomplete elements often lack detectable similarity to
  known TEs, so they're exactly the kind of *bona fide* candidate the filter above would drop.

#### Dependencies

What `slurm/config.sh` actually sets up, and the variable to edit for each if your cluster
differs:

- **RepeatModeler** — run from a shared Apptainer/Singularity container, not installed
  directly: `TETOOLS_CONTAINER_CMD` (`apptainer`) + `TETOOLS_SIF`
  (`dfam-tetools-latest.sif`).
- **BLAST+ / SeqKit / CD-HIT** — loaded as cluster modules in `load_filtering_environment`
  (`ncbi-blast/2.12.0`, `seqkit/2.9.0`, `cd-hit`).
- **MCHelper** — run inside a micromamba env (`ENV_MCHELPER`, `mchelper` by default), pointed to
  by `MCHELPER_PY`. Download it from its repository: <https://github.com/GonzalezLab/MCHelper>,
  then set `MCHELPER_PY` to your own install. The micromamba paths, `MCHELPER_PY`,
  `REPBASE_NT_FASTA` and `TE_PEP_FASTA` are set in the project's `pipeline_config.sh` (see the
  [main README](../README.md#getting-the-pipeline-into-a-project)).
- **BUSCO lineage** — `BUSCO_LINEAGE` in `pipeline_config.sh`, a path to a downloaded lineage
  dataset; species-specific, must be set for your own species.

Run `slurm/check_dependencies.sh` first — it checks every one of these against your edited
`config.sh` before the pipeline itself runs anything.

#### Running it on a SLURM cluster

The SLURM version of the pipeline lives in [`slurm/`](slurm/) — cloning the separate
[TE Candidate Discovery](https://github.com/adriludwig/TE-candidate-pipeline) repository is not
required. To run it on your own project:

1. Edit [`slurm/config.sh`](slurm/config.sh): the cluster module / micromamba / container
   commands, the BUSCO lineage for your species, the curated TE library, the RepBase and
   TE-peptide paths, and the `#SBATCH` resources in each `Step*.slurm.sh` header.
2. Launch the **full pipeline** — this works from any directory, because
   `run_pipeline_slurm.sh` `cd`s into `02_TE_annotation_Clean_Up/` itself before submitting:
   ```bash
   bash 02_TE_annotation_Clean_Up/slurm/run_pipeline_slurm.sh
   ```

`run_pipeline_slurm.sh` runs the dependency check, builds the genome table from `genomes/*.fna`,
then submits Steps 1–4 in order with `sbatch --wait`.

**To re-run one step**, submit it **from inside `02_TE_annotation_Clean_Up/`** — the
individual `StepN_*.slurm.sh` scripts resolve `slurm/config.sh` and their `logs/` folder
from the directory you run `sbatch` in, not from the script's own location:
```bash
cd 02_TE_annotation_Clean_Up
sbatch slurm/Step2_run_TE_candidates_original_headers.slurm.sh
```
Delete that step's output folder first (see [`slurm/README.md`](slurm/README.md) → *Resume
Behaviour*), otherwise the step sees the existing output and skips.

### 3 - Manual Validation of the candidate TEs

This step is essential: the TE Candidate Discovery pipeline's output is a **candidate library**,
and a substantial fraction of those candidates are still expected to be false positives. None
should be added to the curated TE library without going through this manual validation first.
For example, the *T. gondii* genome is devoid of any known TE, even though the pipeline
generates plenty of candidates for it.

Recommended validation steps include:

1. **Inspect TE structure.** Use [TE-Aid](https://github.com/clemgoub/TE-Aid/) (see below) to visualise the
   candidate, inspect its genomic distribution, and evaluate terminal repeats (TIRs or LTRs),
   target site duplications (TSDs), and coding potential. If the TE consensus is still not
   complete, the pipeline can be rerun, excluding step 1 (RepeatModeler) and increasing the
   size of flanking regions or the number of iterative extensions for MCHelper. A manual
   consensus can also be retrieved by hand for specific candidates, following
   [Goubert et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8969392/).
2. **Search for conserved protein domains.** Analyse predicted ORFs using NCBI CD-Search to identify domains characteristic of transposable
   elements (e.g. reverse transcriptase, integrase, RNase H, transposase, endonuclease).
3. **Compare with known TE families.** Search the candidate against RepBase, Dfam, or other TE
   databases to determine whether it represents a highly divergent member of a known family or
   a potentially novel family. An easy search is the online search against the Repbase (https://www.girinst.org/censor/index.php)
4. **Verify TE classification.** Although MCHelper provides an initial structural
   classification, this assignment should always be manually confirmed — it can assign an
   incorrect order or superfamily. Structural features, conserved domains, and sequence
   similarity should all be considered before assigning a final classification.

Only candidates confirmed as bona fide TEs should be added to the curated TE library.


#### TE-aid

One main approach to visualize the candidate TEs is running
[TE-Aid](https://github.com/clemgoub/TE-Aid/) — one composite figure per candidate
(genome-coverage plot, self dot-plot, ORFs, hits to known TE proteins), used to decide
whether a consensus is a *bona fide* TE and whether it is full-length.

Script: [`scripts/run_te_aid.slurm.sh`](scripts/run_te_aid.slurm.sh). TE-Aid takes **one
consensus per run**, so the wrapper splits the candidate FASTA and, for each candidate, runs
TE-Aid against the genome it was recovered from (candidate headers are
`<genome_id>|<candidate_id>`, so the wrapper knows which assembly in `genomes/` to pass to
`-g`).

**Dependencies:** [TE-Aid](https://github.com/clemgoub/TE-Aid/) itself (`TE_AID` in
`pipeline_config.sh` points to your own clone/install), R, NCBI BLAST+ and hmmer — micromamba env
`te_annot` (`ENV_TE_ANNOT`).

Adapt at the top of the script:

- `INPUT_FASTA` / `OUTDIR` — the candidate FASTA to inspect, set as `INPUT_FASTA` + `OUTDIR`
  pairs at the top of the script. Paths are relative to the folder you submit from (the
  project root), and the candidate-pipeline outputs live under `02_TE_annotation_Clean_Up/`,
  so the defaults carry that prefix — adjust it if yours are elsewhere. **The default input
  is `…/final_TE_candidates_after_MCHelper/final_potential_new_TEs.fa`** (the filtered set).
  To inspect the unfiltered non-redundant set instead
  (`all_MCHelper_curated_sequences_NR.cdhit80.fa` — shows the candidates the RepBase/peptide
  filter dropped, e.g. solo-LTR-only families), **comment the default pair and uncomment that
  one**. Each pair writes its own `OUTDIR`, so you can run the script once per file without
  overwriting anything.
- `GENOMES` / `GENOME_SUFFIX` — the step-01 assemblies folder (default `genomes` / `.fna`).
- the `### ADAPT ###` block — put R, NCBI BLAST+ and hmmer on the PATH (the `te_annot` env).
- the `#SBATCH` lines.

**Run it from your project folder** (the one that holds `genomes/`):

```bash
sbatch 02_TE_annotation_Clean_Up/scripts/run_te_aid.slurm.sh
```

**Output:** one `<OUTDIR>/<candidate>/` folder per candidate with its TE-Aid figure, plus
`<OUTDIR>/te_aid_run_summary.tsv` (which genome each candidate ran against, and whether it
was found).

Observations from the first fungi genomes studied: 

- Number of candidate TEs can vary widely between species 
- Some TEs are present in only one copy that is not common to other organisms; one explanation might be that the TEs are inserted into Starships that are represented by single copies. Another possibility is that some fungi have mechanism to control TEs called RIP (Repeat-Induced Silencing) that leads to the silencing of TEs in the genome by introducing mutations in the TE sequences. 
- We can find solo LTRs (that are remnants of complete LTR elements where the internal region is missing) without finding a complete counterpart which is uncommon. Thus, we might want to inspect `all_MCHelper_curated_sequences_NR.cdhit80.fa` instead of `final_potential_new_TEs.fa`. 
- Some TEs that have conserved protein domains do not appear in the TE-aid as having similarity to any TE in the protein level; it could be that the peptides are not conserved enough to be detected by the tools used for the TE-aid construction. We can consider update the TE peptides library that TE-aid uses with the bona fine TEs that we are finding in our curation. 
- Interestingly, among the TE consensus predicted for *Pyricularia*, we can see several chimaeras from TE/non-TE sequences and from different TEs due to segmental duplication. Thus, the number of potential candidate TEs increased because cd-hit was not able to resolve. In those cases, we'll need to fish the TEs from the chimaera sequences.  

#### Output: 
`curate_TE_library.fa` - a FASTA library of curated TEs representing the TE repertoire of the species.

### 4 - Complementary TE Screen pipeline

Pipeline available here: [Complementary TE Screen](https://github.com/adriludwig/Complementary-TE-Screen-pipeline)

An aditional pipeline, `TE_Complementary_Screen`, that searches for TE-like regions that may still be missing after a species-specific curated TE library has already been built. It might happen because some TEs are not that repetitive (as we noticed for some species), or they were excluded downsteam in the candidate TE pipeline.

The workflow first masks each genome with the curated TE library. It then runs a second RepeatMasker screen on the already masked genome using a broader DFAM/RepBase-derived library. Hits from the second screen are treated as candidate missing TE fragments, filtered, extracted from the original unmasked genome, dereplicated per genome, dereplicated again globally across all genomes, extended/refined with MCHelper, and merged into a final candidate FASTA for manual inspection.

This pipeline will be time consuming if a large number of genomes (e.g. >15) are analysed. In this case, one can skip this step and only rely on the previous TE library. In this case, it is essential to check for TE-like protein domains in the pangenome results. Alternatively, the broader library could be filtered to contain only TEs from the taxonimically close group of interest, which could reduce the number of sequences to be mathing the genomes. For Fungi species, the Mycomobilome database can be used in this spet.

Note that for either the DFAM curated library and Mycomobilome library, host genes were found within the TE sequences. The issue was reported to DFAM and they will update their database. Thus, not all aditional mathes mean novel TEs. One should always inspect the additional mathes before accepting them as TEs to be included in the curated TE library.

#### Output: 
`final_gap_screen_candidates.cdhit80.fa` - a FASTA library of candidate TEs. These go through
the same [manual validation](#3---manual-validation-of-the-candidate-tes) as the TE Candidate
Discovery pipeline's own output before anything from this file is added to the curated TE
library.

[↑ Back to top](#table-of-contents)

---

## TE annotation

Once we obtain the curated TE library, we use it to annotate the TEs in every genome with
RepeatMasker.

Script: [`scripts/run_repeatmasker.slurm.sh`](scripts/run_repeatmasker.slurm.sh). It runs
RepeatMasker on every `*.fna` assembly and writes one folder per genome, `RM_<genome_id>/`
(`<genome_id>` = the file name without `.fna`; trimmed to the `GCA_…`/`GCF_…` accession if
the name is one). `-pa` is set to `cpus / 4` (rmblast uses 4 cores per unit).

**Dependencies:** RepeatMasker — micromamba env `repeatmasker4.2`. See
[Installing RepeatMasker](#installing-repeatmasker) below for the Dfam setup it needs.

Adapt at the top of the script:

- `GENOMES="genomes"` — the assemblies folder. By default the `genomes/` folder built in
  **step 01**; change it (or give a full path) if your assemblies are elsewhere.
- `LIB="curated_TE_library.fa"` — the curated TE library FASTA, expected in the project
  folder; change it (or give a full path) if it is elsewhere.
- the `### ADAPT ###` block — put `RepeatMasker` on the PATH (micromamba / module).
- the `#SBATCH` lines — cpus / mem / time / partition.

**Run it from your project folder** (the one that holds `genomes/`):

```bash
sbatch 02_TE_annotation_Clean_Up/scripts/run_repeatmasker.slurm.sh
```

#### Installing RepeatMasker

```bash
micromamba install -c conda-forge -c bioconda repeatmasker
```

RepeatMasker 4.x needs a **Dfam library** configured before it will start. The conda package
does **not** bundle it, so out of the box RepeatMasker aborts with:

```
famdb.py … FamDB data directory not found
```

Fix it once:

1. **Get the Dfam data** if you don't have it — download at least the **root partition (`.0`)**
   from <https://www.dfam.org/releases/> → `families/FamDB/`, `gunzip` the `.h5.gz` into a
   directory (partition 0 alone is enough for `-lib` runs; multi-GB).
2. **`perl $CONDA_PREFIX/share/RepeatMasker/configure`** — accept the TRF default, pick
   **RMBlast** (its `bin` dir), then `5`. For the FamDB prompt, give the directory that holds
   `famdb.py` (`…/share/famdb-*/`).
3. **Point famdb at the data.** `famdb.py` looks for the `*.h5` files in `FAMDB_DATA_DIR` (set
   in `…/share/famdb-*/famdb.conf`), else in `Libraries/famdb` next to `famdb.py`. If your
   files are in a differently-named directory (e.g. `Libraries/Dfam`), edit `famdb.conf`:

   ```
   [famdb]
   FAMDB_DATA_DIR = /full/path/to/the/dir/with/dfamNN.*.h5
   ```

4. **Verify** (from outside the RepeatMasker directory):

   ```bash
   famdb.py info          # prints "Dfam <version>"
   RepeatMasker -v        # version + Dfam, no famdb.py error
   ```

If `famdb.py` reports `ModuleNotFoundError: h5py`, `micromamba install -c conda-forge h5py`
and redo step 4.

> [!TIP]
> The **Dfam / TETools container** (`dfam-tetools-latest.sif`) has all of this pre-configured.
> To use it instead, swap the `### ADAPT ###` block for
> `RepeatMasker() { apptainer exec -B <mount> <sif> RepeatMasker "$@"; }` — but note that
> RepeatMasker's internal file copy can fail on Lustre/NFS from inside a container
> (`cp: cannot lseek … File descriptor in bad state`); if so, stage each genome in node-local
> scratch first.

[↑ Back to top](#table-of-contents)

---

## TE composition analysis

After running RepeatMasker across all genomes, the results can be summarized to evaluate total TE abundance, superfamily distributions, and family-level composition per genome.

Scripts: [`src/pangenome_pipeline/step02_te_annotation_clean_up/summarize_te_content.py`](../src/pangenome_pipeline/step02_te_annotation_clean_up/summarize_te_content.py) ·
[`scripts/summarize_te_content.slurm.sh`](scripts/summarize_te_content.slurm.sh).

**Dependencies:** python3 + pandas + matplotlib — micromamba env `synteny`, already installed
if you set it up for [step 01](../docs/step01_genome_selection_qc.md#dependencies).

#### What this script does:
1. **Filtering:** Drops non-TE hits (such as `Simple_repeat` and `Low_complexity`) from each RepeatMasker `.out` file, writing cleaned tables to `<outdir>/filtered_out/<genome_id>.filtered.out`.
2. **Non-redundant base calculation:** Merges overlapping RepeatMasker hits per contig so overlapping fragments are counted exactly once.
3. **Multi-level summaries:** Computes per-genome total TE content (bp and %), superfamily composition, and family-level breakdown, as well as per-contig statistics.
4. **Data Visualization:** Generates publication-ready summary charts:
   - `genome_size_bp.png` & `total_te_bp.png`: Overall genome sizes and TE amounts.
   - `total_te_pct.png`: Total TE percentage per genome.
   - `superfamily_composition.png`: Stacked bar chart showing relative proportions of TE superfamilies.
   - `family_composition.png`: Stacked bar chart of top TE families.
   - `divergence_vs_total_te.png`: Scatter plot evaluating sequence divergence from consensus vs. total TE content.

#### Usage:

Run it **from your project folder** — the one that holds the `RM_<genome_id>/` folders. It is
a light job (a couple of cores, a few GB); the SLURM wrapper has an **`ADAPT`** block at the
top (`#SBATCH` resources, the env line for a `python3` with `pandas` + `matplotlib` - the
`synteny` env, and the `INDIR` / `OUTDIR` / top-N variables).

```bash
sbatch 02_TE_annotation_Clean_Up/scripts/summarize_te_content.slurm.sh
```

**Parameters:**
- `--indir`: Directory containing the RepeatMasker output subdirectories (`RM_<genome_id>/` or `R_<genome_id>/`, each containing one `.out` and one `.tbl` file — genome length is read from the `.tbl`).
- `--outdir`: Directory where summary TSV tables, filtered outputs, and plots will be saved.
- `--top-superfamilies`: Number of top superfamilies to display before folding into "Other" (default: `8`).
- `--top-families`: Number of top families to display before folding into "Other" (default: `15`).

#### Output files:
- `te_content_summary.tsv`: Global TE metrics per genome (total bp, % genome, fragment count, avg length, divergence).
- `superfamily_summary.tsv` & `superfamily_composition_table.tsv`: Superfamily-level breakdown per genome.
- `family_summary.tsv` & `family_composition_table.tsv`: Family-level breakdown per genome.
- `per_contig/<genome_id>.te_per_contig.tsv`: Detailed TE coverage per contig/scaffold.
- `plots/`: Directory containing all generated `.png` figures.

[↑ Back to top](#table-of-contents)

---

## Filtering TEs from gene annotations

Automated gene prediction pipelines can misclassify TE coding sequences (e.g. transposases, reverse transcriptases) as protein-coding genes. In a pangenome analysis, these spurious gene models often appear as "cloud" genes and artificially bias functional enrichment results.

To resolve this, we can identify TEs that overlap with gene models using the script below and filter the GFFs. Deciding which TEs to filter is not a simple task, though, but an be done by looking at the overlap percentage. From a pilote run with *P. oryzae* and *A. fumigatus* genomes, it seems good cut-offs are:
- **`< 20%` (Keep):** Minor overlap — usually just an edge-clipped or degenerate TE fragment brushing the gene model, not evidence that the gene itself is TE-derived. Kept as a normal protein-coding gene.
- **`20–80%` (Review):** Intermediate overlap, could go either way — a real TE-derived model, or a real gene that happens to sit on/near a TE. Written to `flagged_genes.tsv`/`.gff`/`_cds.fna` — **mandatory manual check** before deciding either way (see [Quantifying TE Overlap Across Gene Models](#quantifying-te-overlap-across-gene-models)).
- **`≥ 80%` (Eliminate):** High-confidence TE-derived gene model. Written to `eliminated_genes.tsv`/`.gff`/`_cds.fna`, then actually removed from the gene set via a kill list + AGAT (see [Removing TE-Derived Genes with a Kill List](#removing-te-derived-genes-with-a-kill-list)). Within this band, genes whose dominant TE match is also unusually divergent (`TE_divergence_pct` above the review cutoff — possible domesticated gene or a weaker call than the overlap % alone suggests) are split out into `eliminated_genes_review.*` for a manual check before trusting the automatic call.

---

### Quantifying TE Overlap Across Gene Models

Scripts: [`scripts/te_gene_overlap.slurm.sh`](scripts/te_gene_overlap.slurm.sh) (SLURM
wrapper) · [`src/pangenome_pipeline/step02_te_annotation_clean_up/te_gene_overlap.py`](../src/pangenome_pipeline/step02_te_annotation_clean_up/te_gene_overlap.py) ·
[`src/pangenome_pipeline/step02_te_annotation_clean_up/extract_flagged_genes.py`](../src/pangenome_pipeline/step02_te_annotation_clean_up/extract_flagged_genes.py) ·
[`src/pangenome_pipeline/step02_te_annotation_clean_up/extract_cds_sequences.py`](../src/pangenome_pipeline/step02_te_annotation_clean_up/extract_cds_sequences.py) (what the wrapper
calls, in that order).

**Dependencies:** python3 + pandas + matplotlib — micromamba env `synteny`, already installed
if you set it up for [step 01](../docs/step01_genome_selection_qc.md#dependencies).

One job does the whole review, writing everything into a single folder:

1. **`te_gene_overlap.py`** intersects RepeatMasker annotations with the GFF gene models
   across all genomes. For protein-coding genes, it calculates the non-redundant TE overlap
   strictly across the union of **CDS** intervals. Every other per-gene fact you'd want for
   review — description, coordinates, dominant TE, `TE_divergence_pct`, overlap position — is
   already a column in its output table, so there's normally no need to go back to the GFF at
   all.
2. **`extract_flagged_genes.py`** + **`extract_cds_sequences.py`**, run three times for three
   non-overlapping bands of `gene_te_overlap.tsv` — writing, for each, a `<prefix>.tsv`
   (same columns, just that subset), a per-genome `<genome_id>.<prefix>.gff` (the genes'
   original GFF records), and a `<prefix>_cds.fna` (spliced CDS FASTA):
   - **`flagged_genes`**: `MIN_PCT` ≤ overlap < `MAX_PCT` (the "inspect" band) — mandatory
     check before eliminating any of these.
   - **`eliminated_genes`**: overlap ≥ `MAX_PCT` — every gene the overlap call alone would
     eliminate.
   - **`eliminated_genes_review`**: overlap ≥ `MAX_PCT` **and** `TE_divergence_pct` ≥
     `REVIEW_MIN_DIV` — eliminated, but on a divergent and/or fragmented TE match (see
     `TE_divergence_pct` below); a high overlap % with a low-confidence match can be a
     domesticated gene rather than a true TE, so these are worth a manual look before
     trusting the elimination.

Run it **from your project folder**. The wrapper has an **`ADAPT`** block at the top
(`#SBATCH` resources, the env line for a `python3` with `pandas` + `matplotlib` - the `synteny`
env, and these variables):

| variable | meaning |
|---|---|
| `TE_INDIR` | folder with the `RM_<genome_id>/` RepeatMasker output dirs, one `.out` each (default `.`) |
| `GFF_INDIR` | folder with one GFF3 per genome (default `gff_originals`) |
| `GENOMES` | folder with one assembly FASTA per genome — only needed for the CDS-extraction steps (default `genomes`) |
| `OUTDIR` | where every output lands (default `te_annotation_integration`) |
| `GENOME_MAP` | optional 2-column TSV (`genome_id<TAB>gff_path`) to force a genome-to-GFF match when the automatic matching (below) gets one wrong |
| `MIN_PCT` / `MAX_PCT` | the "review" band's coverage bounds: `TE_overlap_pct` (% of the gene's CDS a TE covers) between these two lands a gene in `flagged_genes.*`; `MAX_PCT` also marks where "eliminated" starts (overlap ≥ `MAX_PCT`). Default `20`/`80` |
| `REVIEW_MIN_DIV` | minimum `TE_divergence_pct` for an eliminated gene to also land in `eliminated_genes_review.*`. Default `30` |
| `EXTRA_GENES` | optional 2-column TSV (`genome_label<TAB>gene_id`), genes to add to the "review" band's CDS extraction on top of `MIN_PCT`/`MAX_PCT` |

```bash
sbatch 02_TE_annotation_Clean_Up/scripts/te_gene_overlap.slurm.sh
```

Each `RM_<genome_id>/` is matched to its GFF by name (exact / prefix / accession,
case-insensitive) and, failing that, by which GFF's seqids overlap the RepeatMasker contigs
— so mixed naming (accession folders + `AF293` / `A1160` / … folders) is handled
automatically; `GENOME_MAP` above is the manual override for when that automatic match is
wrong.

**Only protein-coding genes and pseudogenes are considered.** RNA-coding genes (tRNA, rRNA,
ncRNA, …) are not.

**The span used:** `protein_coding` = the union of its CDS intervals. `pseudogene` = CDS
union if it has one, else exon union, else the gene's own genomic span (pseudogenes
commonly have no CDS at all).

Every excluded gene is counted and reported at the end of the run:
```
WARNING: N gene(s) had a gene_biotype that isn't 'protein_coding' or 'pseudogene' …
```
A small number here (a handful of real tRNA/rRNA genes) is normal — ignore it. A huge number
means real protein-coding genes are being skipped by mistake — go check that genome's GFF
(how the category is read: an explicit `gene_biotype`/`biotype`/`gene_type` attribute if the
GFF has one, else the feature type in column 3 — see `infer_biotype()` in
`te_gene_overlap.py` for the exact rules).

**Key outputs, all in `OUTDIR` (`te_annotation_integration/` by default):**
- `gene_te_overlap.tsv`: One row per gene **with any TE overlap** (`TE_overlap_bp > 0`) — the
  ~99.7% of genes at 0% are left out of the file (they'd otherwise be nearly all of it) but
  still count toward `n_genes` and every threshold in `protein_coding_threshold_summary.tsv`
  / `pseudogene_threshold_summary.tsv`. Per-gene metrics: `gene_description` (the gene's own
  `description`/`product`/`Note`/`Name` attribute if it has one, else the same from its first
  transcript, whichever the source GFF actually set), `TE_overlap_bp`, `TE_overlap_pct`,
  `TE_divergence_pct` (RepeatMasker's own %-divergence-from-consensus of the TE fragment that
  explains the most of this gene's overlap — a high overlap % on a divergent/fragmented match
  is weaker evidence of a real TE than the same overlap % on a clean, low-divergence one),
  `dominant_TE_class`, `dominant_TE_name`, spatial overlap category `5_prime`, `3_prime`,
  `internal`, or `mixed` (protein_coding genes only).
- `flagged_genes.tsv` / `eliminated_genes.tsv` / `eliminated_genes_review.tsv` (+ matching
  `<genome_id>.<prefix>.gff`): the three bands above — same columns as `gene_te_overlap.tsv`,
  just that subset, plus the matching raw GFF records per genome.
- `multi_transcript_te_overlap.tsv`: Isoform-specific breakdown for multi-transcript protein_coding genes.
- `protein_coding_threshold_summary.tsv` / `pseudogene_threshold_summary.tsv`: Gene count breakdown across candidate overlap cutoffs per genome, one file per category.
- `plots/`: Distribution histograms and threshold sensitivity charts (`protein_coding_*` and `pseudogene_te_overlap_hist.png`).
- `flagged_genes_cds.fna` / `eliminated_genes_cds.fna` / `eliminated_genes_review_cds.fna`: spliced CDS sequence of each band's genes (one record per transcript) — run BLASTn against the TE library or BLASTx/BLASTp against NCBI nr/SwissProt to check whether a flagged model is a real protein or a TE. The header carries `genome|gene_id|transcript`, coordinates, `TE_overlap_pct`, `TE_divergence_pct`, and the dominant TE — cross-reference `gene_id` back into `gene_te_overlap.tsv` for everything else (description, biotype, overlap position, …).
- If a requested gene ends up with no sequence (GFF/FASTA not found, gene or contig missing, no CDS), it's listed in a `WARNING: N requested gene(s) have NO sequence …` block at the end of the run log — check there before assuming the FASTA is complete.

---

### Removing TE-Derived Genes with a Kill List

Scripts: [`scripts/build_te_kill_lists.slurm.sh`](scripts/build_te_kill_lists.slurm.sh) (SLURM
wrapper) · [`src/pangenome_pipeline/step02_te_annotation_clean_up/build_te_kill_lists.py`](../src/pangenome_pipeline/step02_te_annotation_clean_up/build_te_kill_lists.py) (what it calls).

**Dependencies:** python3 + pandas — micromamba env `synteny`, already installed if you set it
up for [step 01](../docs/step01_genome_selection_qc.md#dependencies) — and AGAT
(`agat_sp_filter_feature_from_kill_list.pl`), hardcoded here with an env called `bioinfo`
(adapt to your own).

One job does both steps:

1. Generates a kill list of gene IDs from each genome's own TE-overlap table
   (`gene_te_overlap.tsv`, `TE_overlap_pct` ≥ `MIN_PCT`, default `80`).
2. Runs AGAT (`agat_sp_filter_feature_from_kill_list.pl`) to remove those genes from every GFF
   in `GFF_INDIR`, writing cleaned copies to `GFF_OUTDIR`.

The wrapper has **two** `### ADAPT ###` blocks to edit — one per dependency above — plus these
variables at the top:

| variable | meaning |
|---|---|
| `GENE_TABLE` | path to `gene_te_overlap.tsv` (default `te_annotation_integration/gene_te_overlap.tsv`) |
| `KILL_OUTDIR` | where the kill list is written (default `kill_lists`) |
| `MIN_PCT` | minimum `TE_overlap_pct` to kill a `protein_coding` gene (default `80`) |
| `EXTRA_GENES` | optional 2-column TSV (`genome_label<TAB>gene_id`), extra genes to add to the kill list on top of `MIN_PCT` |
| `GFF_INDIR` | which GFFs to filter — **edit this one:** `gff_originals` (default, the main case — clean before Liftoff) or `liftoff_polished` (if Liftoff already ran; same kill list works, Liftoff preserves the source gene IDs) |
| `GFF_OUTDIR` | where the filtered GFFs are written (default `gff_cleaned`) |

```bash
sbatch 02_TE_annotation_Clean_Up/scripts/build_te_kill_lists.slurm.sh
```

---

## Generating a TE-only GFF3

> 🚧 **To be developed.**

RepeatMasker already runs with `-gff`, so a raw `.out.gff` exists per genome, but it's not
usable as-is: it needs to be filtered (drop `Simple_repeat`/`Low_complexity` hits) and
standardized into proper GFF3. Approach still to be decided.

A standalone TE GFF3 would also be useful on its own for VEuPathDB: loaded as its own track in
JBrowse, separate from the gene models.

[↑ Back to top](#table-of-contents)

---

## Identifying Starships in fungal genomes

> 🚧 **To be developed (Adriana).** 

Starships are giant mobile elements of the Pezizomycotina — large insertions, typically tens
to hundreds of kb (the largest exceed 700 kb) — mobilised by a **captain** tyrosine recombinase
(a gene with a **DUF3435** domain at the 5′ end of the element)
([Gluck-Thaler et al. 2022](https://doi.org/10.1093/molbev/msac109)). Unlike the autonomous
TEs curated above, a Starship carries a **cargo of ordinary host genes** — often many, and
frequently adaptive ones (metal and formaldehyde tolerance, secondary metabolism, virulence
factors).

**Why Starships are *not* filtered out.** The cargo genes are real protein-coding genes with
real cellular functions. Removing Starships from the pangene matrix would discard genuine
accessory content and misrepresent the pangenome. Starships are therefore **kept**. This is
the opposite of [*Filtering TEs from gene annotations*](#filtering-tes-from-gene-annotations)
above, where TE-derived gene models (transposases, reverse transcriptases, …) are removed
precisely because they add no host function.

**What we do need.** Identify the Starships so their signature in the pangenome can be read
for what it is — a block of co-located, co-inherited cargo genes moving as a unit — rather
than as many independent gene gain/loss events. At the end of the pipeline, once the
pangenome is built, we want to:

- delimit Starship insertions in each assembly and list their cargo genes;
- tag those genes in the final pangene table (e.g. `starship_cargo` + element ID), so the
  accessory / cloud fraction can be split into mobile-element-driven vs. other;
- report Starship presence/absence across strains alongside the pangene results.

**Intended approach** (to be finalised): check and adapt
[`starfish`](https://github.com/egluckthaler/starfish) (Gluck-Thaler & Vogan 2024) to find and visualize them.

[↑ Back to top](#table-of-contents) &nbsp;|&nbsp; [Back to top level](../README.md)


