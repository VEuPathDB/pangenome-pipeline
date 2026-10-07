# Pangenome Analysis Pipeline Documentation

This repository documents the **pangenome analysis workflow**. The pipeline consists of several steps, each with its own `README.md` and a `scripts/` folder containing the associated codes.

**Internal document.** This is primarily an internal working document: the scripts are mostly
SLURM scripts, including code already adapted for our shared containers/installations. They
will eventually need editing to add paths to environments that aren't shared. They should also
be adapted/reviewed before being made public, as they contain comments that make sense to our
working group but not to outside readers.


Most pangenome tools call pangenes by
**sequence-similarity clustering**: throw every protein from every genome into one pool,
cluster by similarity, and each cluster becomes one pangene. One of the main issues with the clustering method is that it cannot separate paralogs properly while GET_PANGENES, because it uses the synteny information, paralogs will be considered independent pangenes.

Here we chose to use [GET_PANGENES](https://github.com/Ensembl/plant-scripts), which calls pangenes by **similarity** and **synteny**: it performs whole-genome alignments (WGA) for every assembly to define clusters of collinear, orthologous gene features from each genome's GFF annotation. Such clusters define the pangenes that make up the pangenome
([Contreras-Moreira et al. 2023](https://link.springer.com/article/10.1186/s13059-023-03071-z)).


Both methods are only as good as the gene annotation they're given. Real genes might be missing from the annotation, creating a false absence in the pangenome matrix, and false gene calls might inflate the accessory genome.

So before running GET_PANGENES, this pipeline was developed to get consistent gene annotation across all genome assemblies:

| Step | |
|---|---|
| [01 – Genome Selection & QC](01_Genome_Selection_QC/README.md) | Selects and quality-checks the genome set itself, and normalizes every assembly/GFF into one consistent naming convention every later step relies on. Profiles every contig's size and whether the small ones carry genes, to see how contiguous each assembly actually is and catch species-specific peculiarities (e.g. Tgondii's fragmented mitochondrial genome). Screens for and removes organelle (mitochondrial/apicoplast) and redundant/duplicate contigs avoiding inflated gene counts from genome redundancy, and false PAV from inconsistent organelle-genome inclusion across the panel.|
| [02 – TE Annotation & Clean Up](02_TE_annotation_Clean_Up/README.md) | Builds a species-specific curated TE library and annotates every genome with it, giving a TE-composition picture per genome (total content, superfamily/family breakdown). Scores how much of each annotated gene's CDS overlaps a TE and filters out genes that are probably TEs. Keeping TEs annotated as genes inflates the cloud (rare/strain-specific) gene category, since TEs are expected to insert at different locations in different strains. **Starships** (giant mobile elements in fungi that carry real host genes as cargo) are deliberately kept, not removed — their cargo genes are genuine; identifying the Starships is a step we need to do later when we have the pangenome results. |
| [03 – *de novo* Genome Annotation](03_denovo_Genome_Annotation/README.md) | **Optional** — only for genomes with no gene annotation at all. For fungal genomes, runs **Funannotate**, a standard *de novo* annotator. Without this step, the assembly could still get by on genes transferred from other strains (step 05's Liftoff) — but any gene truly unique to this strain, with no counterpart in any other genome's annotation, would never be captured that way. |
| [04 – InterProScan](04_InterProScan/README.md) | Some genomes may already carry InterPro domain annotation in their description, but for every genome to be on the same InterProScan version, this step re-runs it across the panel — the presence of an InterPro domain is used as evidence in step 05. |
| [05 – Cross-genome annotation integration](05_Cross-genome_annotation_integration/README.md) | Different genomes usually arrive annotated by different groups/tools, with subtly different gene-calling tools/conventions. To get consistent annotation across all genomes, this step re-annotates every genome. Our "all-vs-all" strategy has several steps starting with the LiftOff of all strains against all strains. To build the final annotation, we have 2 approaches. **Approach 1: filling each locus by a fixed priority order of sources**  1) pick a trusted reference genome, transfer its genes onto every other genome with Liftoff (keeping only high-confidence transfers — valid gene model, high enough coverage and sequence identity); 2) fill gaps from that genome's own original annotation; 3) then fill whatever's still missing from every other genome's annotation too (also via Liftoff, high-confidence transfers); 4) then fill whatever's still missing from every other genome's annotation (low-confidence). After obtaining the final annotations, one important step is to evaluate the necessity of  some filtering. Overlapping genes arise from the pipeline and need to be filtered, keeping only the more confident gene model. Also, some very low-quality transfers get added to the GFFs in the last steps — some may remain, but others should be filtered out entirely (e.g. a gene transfer with <50% coverage). **Two alternative approaches to this filtering are documented too** — resolving and filtering before `get_pangenes` runs, or adding specific tags to gene description so we can address the low quality genes after running GET_PANGENES. **Approach 2: filling each locus by whichever transferred gene model is individually best-supported regardless of source**. SADIK TO WRITE THE SUMMARY HERE.<br><br>**Caveats:**<br>- a transferred gene model having a valid ORF doesn't mean it's expressed or functional, since the pipeline only checks the coding sequence itself — a mutation in a regulatory region (e.g. promoter) that silences the gene would go undetected;<br>- gene models with invalid ORFs, called "pseudogene" in the original GFFs ot the ones kept in our final GFFs (those transferred with invalid ORFs) could be genuine pseudogenes or may actually be functional genes whose CDS merely carries a sequencing error. Keeping them in get_pangene analysis minimizes the chance of false PAV call. |
| [06 – GET_PANGENES](06_GET_PANGENES/README.md) | Links each genome's already-ready GFF and FASTA into one directory, paired up by filename the way GET_PANGENES expects, then builds the actual pangenome — the similarity + synteny gene calling described above — from the cleaned-up, cross-genome-consistent annotations steps 01–05 produced. |
| [07 – Post-GET_PANGENES adjusting](07_Post_GET_PANGENES_adjusting/README.md) | GET_PANGENES sometimes merges two adjacent genes into one pangene call. This step is to understand them and separate to distinct pangenes when necessary.   |



## Steps

- **[01 - Genome Selection & QC](01_Genome_Selection_QC/README.md)** 
- **[02 - TE Annotation & Clean Up](02_TE_annotation_Clean_Up/README.md)**
- **[03 - *de novo* Genome Annotation ](03_denovo_Genome_Annotation/README.md)** 
- **[04 - InterProScan](04_InterProScan/README.md)** 
- **[05 - Cross-genome annotation integration](05_Cross-genome_annotation_integration/README.md)** 
- **[06 - GET_PANGENES](06_GET_PANGENES/README.md)** 
- **[07 - Post-GET_PANGENES_adjusting](07_Post_GET_PANGENES_adjusting/README.md)** 
- **[08 - Pangenome Visualisation](08_Pangenome_Visualisation/README.md)** 


---

## Pangenome analysis reports

Results of this pipeline for each species analysed so far. These are curated reports (the full
technical write-up plus the figures/tables it references) - the full working data behind them
(raw assemblies, intermediate GFFs, InterProScan/Liftoff output, etc.) is not published here, as
it runs from hundreds of MB to several GB per species.

- **_Aspergillus fumigatus_**
  - [Summary](results/Afum/SUMMARY.md)
  - [Full report](results/Afum/README.md)
- **_Toxoplasma gondii_**
  - [Summary](results/Tgondii/SUMMARY.md)
  - [Full report](results/Tgondii/README.md)

---

## Alternative Get-PanLoci Pipeline

GET_PANGENES calls pangenes from annotated gene models, so its results are bounded by annotation quality — the reason steps 01–05 exist. The [Alternative Get-PanLoci Pipeline](Alternative_Get_PanLoci_pipeline/README.md) takes a different approach for cases where annotation is inconsistent or missing across the panel: it builds a whole-genome alignment with Cactus and uses `halLiftover` to project each **panlocus** (a genomic interval, not necessarily an annotated gene) across every genome, then validates each transfer at the DNA sequence level with `minimap2`, independent of whether the target strain has a gene call there. It runs as its own pipeline rather than as a step in the sequence above.

---

## Getting the pipeline into a project

Clone or download this repository into the root folder of the species project you're working
on, so its step folders sit directly under it, alongside your own `genomes/`,
`gff_originals/`, and `genome_order.txt`:

```bash
git clone https://github.com/adriludwig/pangenome-pipeline.git <project>/
```

```
<project>/
├── 01_Genome_Selection_QC/
│   └── scripts/
├── 02_TE_annotation_Clean_Up/
│   └── scripts/
├── 03_denovo_Genome_Annotation/
│   └── scripts/
├── 04_InterProScan/
│   └── scripts/
├── 05_Cross-genome_annotation_integration/
│   └── scripts/
├── 06_GET_PANGENES/
│   └── scripts/
├── 07_Post_GET_PANGENES_adjusting/
│   └── scripts/
├── 08_Pangenome_Visualisation/
│   └── scripts/
├── lib/pangenome_common/       Python helpers shared by the step scripts (keep it with them)
├── species_specific/           one-off scripts for a single species (Afum/, Tgondii/)
├── genome_order.txt            the genome list (generated in step 01)
├── genomes/                    <basename>.fna    – assemblies
└── gff_originals/              <basename>.gff3   – each genome's own annotation
```

`genomes/`, `gff_originals/`, `genome_order.txt` are yours to create in step 01. 

**Cluster settings - once per project.** Where micromamba, the tool installs and the databases live
on your cluster is set in one file, `pipeline_config.sh` at the project root. It is not part of the
repository; start from the template:

```bash
cp pipeline_config.example.sh pipeline_config.sh   # then fill in the paths you need
```

Every job script loads it through [`lib/slurm_common.sh`](lib/slurm_common.sh), so submit jobs from
the project folder. A job stops with a message naming any setting it needs that is still empty.
Any setting can also be given for a single job, e.g. `sbatch --export=ALL,LABEL=Toxo7a ...`.

What stays in each script: the `#SBATCH` lines, the `### ADAPT ###` block (cluster module names;
the micromamba environment it activates is named in `pipeline_config.sh`), and the run's own
settings at the top of the script (input folders, thresholds, and species-specific values such as
`SPECIES`, `LABEL` or the organelle contigs, which have no default and must be set).

If a script fails with `Permission denied`, restore the executable bit on all scripts:

```bash
chmod +x $(find . -name "*.sh")
```