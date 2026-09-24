# 01 – Genome Selection & QC

***Purpose:*** Select the genomes for the pangenome analysis; confirm each one's species
identity and check its assembly for organellar, redundant, or otherwise unwanted contigs; and
put every assembly/GFF into a single, consistently named form that every later step relies on.

## Dependencies

| Program | Used by |
|---|---|
| skani | `check_species.slurm.sh` |
| minimap2 | `synteny_vs_ref.slurm.sh`, `contig_redundancy.slurm.sh` |
| samtools | `contig_size_profile.slurm.sh`, `contig_redundancy.slurm.sh`, `check_redundant_genes.slurm.sh`, `screen_organelles.slurm.sh` |
| blastn (BLAST+) | `check_redundant_genes.slurm.sh`, `screen_organelles.slurm.sh` |
| python3 + pandas + matplotlib | `contig_size_profile.slurm.sh`, `synteny_vs_ref.slurm.sh` |

These scripts are hardcoded with the author's own micromamba environment (paths and env name
`synteny`) - edit the `### ADAPT ###` block at the top of each script to point to your own.
Easiest fix: create one env with the same name and every program above installed:

    micromamba create -n synteny -c bioconda -c conda-forge skani minimap2 samtools blast python pandas matplotlib scipy seaborn

`scipy` and `seaborn` are only needed by the plots of [step 08](../08_Pangenome_Visualisation/README.md). If `synteny` already exists, add them with
`micromamba install -n synteny -c conda-forge scipy seaborn`.


## Table of Content

1. [Issues addressed by this step](#1-issues-addressed-by-this-step)
2. [Directory layout](#2-directory-layout)
3. [Obtaining the genomes](#3-obtaining-the-genomes)
4. [Naming convention](#4-naming-convention)
5. [Trim the FASTA headers](#5-trim-the-fasta-headers)
6. [The genome list](#6-the-genome-list-genome_ordertxt)
7. [Species check](#7-species-check)
8. [Contig size profile](#8-contig-size-profile)
9. [Synteny vs the reference](#9-synteny-vs-the-reference)
10. [Organelle screen](#10-organelle-screen)
11. [Contig-redundancy screen](#11-contig-redundancy-screen)
12. [Remove the unwanted contigs](#12-remove-the-unwanted-contigs)

---

## 1. Issues addressed by this step

### 1.1 Species identity

Public assemblies are sometimes mislabelled. One wrong-species genome inflates the accessory pangenome. Confirm
every genome against the reference with whole-genome ANI ([species check](#7-species-check),
section 7) — anything below ~95% ANI is not the same species; drop it or investigate.

### 1.2 Organellar sequences

Some assemblies include non-nuclear genomes — the **mitochondrion**, and depending on the
lineage the **apicoplast** (Apicomplexa). If some genomes carry an organelle with its genes annotated
in the GFF and others don't, those genes become false presence/absence in the pangenome.
Check every genome ([organelle screen](#10-organelle-screen)) and strip the organelle contigs
([section 12](#12-remove-the-unwanted-contigs)) so the whole set is nuclear-only and consistent.


### 1.3 Assembly fragmentation and small contigs

GET_PANGENES (step 06) calls genes from synteny, so a very fragmented assembly gives poor
results — and short contigs are often duplicated pieces of a chromosome (uncollapsed repeats,
haplotig-style overlaps, extra rDNA or organelle copies) whose genes, if redundant, get
counted as spurious accessory genes and inflate the accessory pangenome.

Start with the [contig size profile](#8-contig-size-profile) (section 8): a cheap first look
at contig count, size, and gene content for every genome. Use it to drop or flag an assembly
that's too broken overall (the [synteny dot-plots](#9-synteny-vs-the-reference) help confirm
this visually) and to shortlist which small contigs are worth checking further. Only then run
the [contig-redundancy screen](#11-contig-redundancy-screen) (section 11), which does the
actual duplicate check: it flags small contigs that align to a longer contig of the same
assembly, you check the genes on them, and the confirmed duplicates go on the drop list
([section 12](#12-remove-the-unwanted-contigs)).

> ⚠️ **WARNING.** Not every small/redundant-looking contig is an assembly artefact — two
> different issues can produce the same pattern and need to be told apart before dropping
> anything:
>
> - **Genuine mini/accessory chromosomes (fungi).** Several fungal species carry real mini /
>   accessory chromosomes (~0.5–3 Mb): *Pyricularia oryzae*
>   ([Langner et al. 2021](https://doi.org/10.1371/journal.pgen.1009386)), *Fusarium oxysporum*
>   ([Yang et al. 2020](https://doi.org/10.1094/PHYTO-03-20-0069-IA);
>   [Zhang et al. 2020](https://doi.org/10.1038/s42003-020-0770-2)), *Zymoseptoria tritici*
>   ([Goodwin et al. 2011](https://doi.org/10.1371/journal.pgen.1002070)), *Fusarium solani*
>   ([Coleman et al. 2009](https://doi.org/10.1371/journal.pgen.1000618)), and others. At least
>   in *P. oryzae* it is described that minichromosomes originate by segmental duplication and
>   rearrangement of core chromosomes ([Langner et al. 2021](https://doi.org/10.1371/journal.pgen.1009386)).
>   We need to explore the literature to understand if we can separate genuine minichromosomes
>   from redundant contigs that are assembly artefacts or decide if eliminating them is the best approach.
> - **Fragmented mitochondrial genome (*T. gondii*).** The *T. gondii* mitochondrion does not
>   assemble as one clean contig — it comes out fragmented across several short contigs/pieces ([Namasivayam et al., 2021](https://genome.cshlp.org/content/31/5/852)).

### 1.4 Contamination

⚠️ **Not yet addressed in this pipeline:** We assume that the steps of identifying and eliminating contamints were already done by the authors of the genomes; however, true contamination (human/bacterial DNA) could still be present in a few of them. e.g. During the investigation of repeat sequences identified by RepeatModeler in one of *P. falciparum* strains, we identified a human TE sequence and realized that 2 contigs were from a contamination origin. The same issue was found in *T. gondii*: TE-candidate curation identified human L1/LINE-1 sequences in multiple small contigs/strains, including the reference genome ME49.

[↑ Back to top](#table-of-content)

---

## 2. Directory layout

This step needs two folders at the project root — `genomes/` and `gff_originals/` — alongside
`genome_order.txt` (section 6); see the
[full project layout](../README.md#getting-the-pipeline-into-a-project) in the main README for
where every step's own folder sits relative to these. Create the two here:

```bash
mkdir -p genomes gff_originals
```

Run every command below **from `<project>/`** — so `genomes/`, `gff_originals/`,
`genome_order.txt` are in the current directory and the scripts are at
`01_Genome_Selection_QC/scripts/`.

Keep the two folders separate and the basenames identical (`genomes/Af293.fna` ↔
`gff_originals/Af293.gff3`). Downstream that is `--genomes-dir genomes` /
`--gff-indir gff_originals` (step 02), the InterProScan inputs (step 04), the Liftoff inputs
(step 05), and so on. Genomes annotated *de novo* in step 03 get their FASTA in `genomes/`
now and their GFF added to `gff_originals/` at the end of step 03.


[↑ Back to top](#table-of-content)

---

## 3. Obtaining the genomes

For every genome: an assembly FASTA and its annotation GFF (unless it will be annotated *de novo* in step 03).

### 3.1 From NCBI

Install the NCBI `datasets` tool:

```bash
micromamba install -c conda-forge ncbi-datasets-cli
```

Pass every accession in one call, then sort the archive into the two folders:

```bash
mkdir -p genomes gff_originals

datasets download genome accession \
    GCA_000002495.2 \
    GCA_004346965.1 \
    GCA_049355315.1 \
    ...                        # one line per accession
    --include genome,gff3

unzip ncbi_dataset.zip
# datasets unpacks into ncbi_dataset/data/<accession>/ - move each pair out, named
# for the accession, then bin the rest
for d in ncbi_dataset/data/*/; do
    acc=$(basename "$d")
    mv "$d"/*_genomic.fna "genomes/${acc}.fna"
    [ -f "$d/genomic.gff" ] && mv "$d/genomic.gff" "gff_originals/${acc}.gff3"
done
rm -r ncbi_dataset ncbi_dataset.zip README.md md5sum.txt
```

`--include genome,gff3` fetches the assembly and its annotation. Result: matching
`genomes/<accession>.fna` / `gff_originals/<accession>.gff3` pairs. Headers still need
trimming (section 5).


### 3.2 From VEuPathDB (FungiDB, PlasmoDB, …) or other sources

No accession-based fetch — download each genome's FASTA and GFF by hand and put them in the
two folders with a **shared basename**: `genomes/<basename>.fna` and
`gff_originals/<basename>.gff3`. The names as downloaded usually don't match (e.g.
`…_TgondiiARI_Genome.fasta` vs `…_TgondiiARI.gff`) — rename both to one short label
(`ARI.fna` / `ARI.gff3`) before continuing. The header-trim check (section 5) catches any mismatch.

[↑ Back to top](#table-of-content)

---

## 4. Naming convention

Each genome is two files with the **same basename** — the assembly and its annotation. The
basename itself doesn't matter.

```
genomes/<basename>.fna
gff_originals/<basename>.gff3
```

**Extensions:** use `.fna` for the assembly and `.gff3` for the annotation. The check and the
step-02 scripts also accept `.fasta`/`.fa` and `.gff`, but keeping the whole set on `.fna` /
`.gff3` avoids surprises later — the download and rename steps already produce that.

[↑ Back to top](#table-of-content)

---

## 5. Trim the FASTA headers

Scripts: [`rename_headers.py`](scripts/rename_headers.py) ·
[`rename_headers_slurm.sh`](scripts/rename_headers_slurm.sh).

`rename_headers.py` needs no other input. For every FASTA it trims each header to its first
token (`>CM000171.1 Aspergillus fumigatus chromosome 1` → `>CM000171.1`) and, if the filename
starts with an NCBI accession, shortens it to just that:

```
GCA_040142875.1_UCR_Afum_08-19-02-30_1.0_genomic.fna  ->  GCA_040142875.1.fna
AfumigatusAf293.fasta                                 ->  AfumigatusAf293.fasta   (unchanged)
```

**`rename_headers_slurm.sh` always edits `genomes/*.fna` in place — that's the only behavior it
has, there's no flag or variable to pick something else at submit time.** If you want copies
written elsewhere instead, you have to open `rename_headers_slurm.sh` and change its
`python3 ... --in-place` line to `--outdir <dir>` yourself (see the argument table below for
what that changes).

Submit it from the project folder, so it finds `01_Genome_Selection_QC/scripts/`, `genomes/`
and `gff_originals/`. It has an ADAPT BEFORE RUNNING block at the top (SBATCH resources, the
module load python line) — set that to your cluster first.

```bash
sbatch 01_Genome_Selection_QC/scripts/rename_headers_slurm.sh
```

It trims every header in genomes/ in place, then cross-checks: every sequence name in each
FASTA also has to appear in that genome's GFF, column 1. Without this, a naming mismatch
(GenBank vs RefSeq IDs, or chr1 vs 1) silently breaks any later step that pulls sequence out of
the genome using its own GFF's coordinates.

Check the job log once it finishes:

- `check-gff: OK` — passed, every ID matched.
- `WARNING: M/N trimmed IDs are not in ...` — some FASTA names don't match the GFF; fix that
  before moving on.

Rename the GFFs so each shares its assembly's basename **before** running this, so the
cross-check has something to compare against.

`rename_headers.py`'s own arguments, for reference (what `rename_headers_slurm.sh` sets for
you, and what you'd change by hand in that file):

| argument | meaning |
|---|---|
| `genomes/*.fna genomes/*.fasta` | the assemblies (any mix of `.fna`/`.fasta`/`.fa`) |
| `--in-place` | edit the files where they are — **this is what `rename_headers_slurm.sh` uses** |
| `--outdir DIR` | write trimmed copies here instead (keeps the originals) |
| `--to-ext .fna` | force every output to this extension (omit it and each output just keeps its input's own extension) |
| `--check-gff gff_originals` | also verify each assembly against its GFF (see above) |

The script itself has no built-in default: called on its own, you must pass exactly one of
`--in-place` / `--outdir`, or it exits with an error
(`ERROR: give --in-place, --outdir, or --check-gff.`). `--check-gff` **by itself**, with
neither of those two, runs the check only — nothing gets written anywhere.

**Output:** `genomes/*.fna` — same files, rewritten in place with trimmed headers. No new
files, nothing added to `gff_originals/`.

[↑ Back to top](#table-of-content)

---

## 6. The genome list (`genome_order.txt`)

List every genome to include, one basename per line, in `genome_order.txt`. Build it now, from
whatever's currently in `genomes/` — every screen starting with section 7 (species check
onward) needs this file to exist before it can run.

```bash
ls genomes/*.fna | sed 's#.*/##; s/\.fna$//' > genome_order.txt
```

Then open the file and put the **reference genome on the first line**. You can reorder the
others too if you want a particular order — it just needs the reference first.

This isn't final. Treat it as a working list: if the [species check](#7-species-check) or the
[contig size profile](#8-contig-size-profile) turns up a genome that shouldn't be in the panel
(wrong species, unusably fragmented), just delete its line here and remove its files from
`genomes/`/`gff_originals/` — nothing downstream depends on the list being locked in at this
point.

**Output:** `genome_order.txt` — one genome basename per line, reference first.

[↑ Back to top](#table-of-content)

---

## 7. Species check

This step, plus [synteny](#9-synteny-vs-the-reference), [organelle](#10-organelle-screen), and
[contig-redundancy](#11-contig-redundancy-screen) later in this file, all draw from the same
toolset: `skani`, `minimap2`, `blast`, `samtools`, and a `python3` with `matplotlib` — e.g. a
micromamba env:

```bash
micromamba create -n synteny -c bioconda -c conda-forge skani minimap2 blast samtools matplotlib scipy seaborn
```

Each `.slurm.sh` has an **`ADAPT BEFORE RUNNING`** block at the top — set the `#SBATCH`
resources and the tool-environment line (micromamba / `module load` / `conda activate`) to
your cluster before submitting. Where a screen has thresholds worth reviewing per species
they are named variables in that block (or `--flags` on the underlying `.py`, listed by
`--help`); the sections below call out which ones actually need attention — the rest have
species-agnostic defaults and can be left alone.

Script: [`check_species.slurm.sh`](scripts/check_species.slurm.sh).

`skani` computes whole-genome ANI of every genome against the reference (line 1 of
`genome_order.txt`).

```bash
sbatch 01_Genome_Selection_QC/scripts/check_species.slurm.sh
```

**Output:** `species_check/species_check.tsv` — one row per genome:

- `ANI_to_<ref>` — the % identity itself. `0.00` specifically means skani couldn't align the
  genome to the reference at all (too divergent), not that it measured 0% identity.
- `aligned_fraction` — how much of the genome skani could actually align to the reference
  (0–1), regardless of how identical those aligned parts are. Informational only — the verdict
  below doesn't use it — but a `same_species` genome with a low aligned fraction is worth a
  look anyway (could be a very incomplete/fragmented assembly).
- `verdict`:
  - `same_species` — ANI ≥ `THRESHOLD` (default **95%**, edit at the top of the script for a
    different genus). Typically ≥ 99% ANI; cryptic siblings can fall to ~90–93%.
  - `CHECK` — below `THRESHOLD`. Either a different species (drop it), or, if the ANI is only
    just under the threshold, a very divergent strain worth a closer look.

[↑ Back to top](#table-of-content)

---

## 8. Contig size profile

Before anything that decides to keep or discard a contig — organelle screening, redundancy
screening, synteny comparisons — get a picture of the raw material first: how many contigs
each assembly has, how big they are, and whether the smallest ones actually carry an
annotated gene. 

Scripts: [`scripts/contig_size_profile.slurm.sh`](scripts/contig_size_profile.slurm.sh)  · [`scripts/contig_size_profile.py`](scripts/contig_size_profile.py).

```bash
sbatch 01_Genome_Selection_QC/scripts/contig_size_profile.slurm.sh
```

For every genome: contig lengths from its `.fai` (built automatically via `samtools faidx` if
missing), plus — if `GFF_INDIR` is set (default: `gff_originals`) — a gene count per contig,
both every gene-level feature (`n_genes`) and, separately, just `protein_coding`/`pseudogene`
(`n_pc_pseudo`) — the biotypes GET_PANGENES actually uses downstream, so a lone tRNA sitting
on an otherwise-empty contig doesn't make that contig look load-bearing for the pangenome.
`GFF_INDIR` and `MIN_LEN` (10,000 bp — the length under which a contig counts as "small") are
both variables at the top of `contig_size_profile.slurm.sh`, edit them there.

**Output,** all in `contig_size_profile/`:

- `contig_sizes.tsv` — every contig, with its length and gene counts.
- `contig_size_summary.tsv` — one row per genome: total contigs, total bp, and how many of
  each fall under `MIN_LEN` (default 10 kb).
- `small_contig_genes.tsv` — one row per gene sitting on a contig under `MIN_LEN`, with its
  description and coordinates. Read this before dropping anything: mitochondrial/apicoplast
  genomes are naturally this small, so a description like "cytochrome c oxidase" or a cluster
  of rRNA/tRNA genes on one contig is the organelle, not contamination. An isolated,
  ordinary-looking nuclear gene on a tiny contig is the case that actually needs a closer look.
- `plots/` — four PNGs: contig count by size bin (kb/Mb, no log scale), small-contig count per
  genome, and gene counts (all genes, then just protein_coding/pseudogene) by size bin.

This is a size/contamination triage, not the redundancy screen
([section 11](#11-contig-redundancy-screen)) — a contig can be short and unique (this
section) or short and a near-duplicate of a longer contig (that one); check both before
deciding a contig's fate.

[↑ Back to top](#table-of-content)

---

## 9. Synteny vs the reference

Scripts: [`synteny_vs_ref.slurm.sh`](scripts/synteny_vs_ref.slurm.sh) ·
[`plot_dotplot.py`](scripts/plot_dotplot.py).

`minimap2` aligns the reference (line 1 of `genome_order.txt`) to every other genome, then one
dot-plot per genome — reference on X, query contigs on Y ordered by where they land on the
reference. Blue = same strand, red = inversion. Edit `PRESET` at the top of
`synteny_vs_ref.slurm.sh` (`asm5`/`asm10`/`asm20`) for more divergent species.

```bash
sbatch 01_Genome_Selection_QC/scripts/synteny_vs_ref.slurm.sh
```

**Output:** `synteny/<label>_vs_<ref>.png` (the dot-plot) and `synteny/<label>_vs_<ref>.paf`
(the raw alignment), one pair per non-reference genome.

[↑ Back to top](#table-of-content)

---

## 10. Organelle screen

Scripts: [`screen_organelles.slurm.sh`](scripts/screen_organelles.slurm.sh) ·
[`screen_organelles.py`](scripts/screen_organelles.py).

Tells you, for each genome, whether the mitochondrion / apicoplast / etc. is present as its
own contig, only as a fragment inserted in a nuclear sequence, or missing.

**Step 1 — identify the organelle contig in one genome you know.** Pick any one genome where
you can tell which sequence is the organelle — from the assembly description, or from the
contig names (VEuPathDB, for instance, names it `mito_...` / `apicoplast_...`; on NCBI the
mitochondrion is flagged in the assembly's sequence report; failing that it's usually the
short ~20–100 kb sequence). You just need its **genome basename** and the **contig name(s)**,
e.g. genome `AF293`, contig `mito_A_fumigatus_Af293`.

**Step 2 — put those in the script and run it.** Edit the two variables at the top of
`screen_organelles.slurm.sh`:

```bash
SOURCE_GENOME="AF293"
ORGANELLE_CONTIGS="mito_A_fumigatus_Af293"      # space-separated if more than one
```

```bash
sbatch 01_Genome_Selection_QC/scripts/screen_organelles.slurm.sh
```

The script pulls those contig(s) out of `SOURCE_GENOME`, `blastn`s them against every genome
in `genome_order.txt`, and `screen_organelles.py` writes **one row per (genome, organelle,
contig) that hits** — so a genome that has the organelle spread over, or *duplicated* across,
several contigs shows up as several rows.

`contig_class` (per contig):

| | |
|---|---|
| `full` | contig ≈ one organelle copy |
| `full_in_scaffold` | whole organelle, but embedded in a large contig |
| `fragment` | partial hit (NUMT / NUPT, or a broken piece) |

`genome_call` (per genome, from all its contigs):

| | |
|---|---|
| `own_contig` | exactly one `full` contig |
| `MULTIPLE_COPIES` | more than one `full` contig — usually an assembly artefact, check before keeping any |
| `split_across_contigs` | no single full contig, but the pieces together cover the organelle |
| `in_scaffold` | a `full_in_scaffold` and no `full` |
| `fragments_only` | only fragments |
| `absent` | no hit |

**Output:** `organelle_screen/organelle_screen.tsv` and a per-organelle summary in the job log.
Use it to build the drop list in [section 12](#12-remove-the-unwanted-contigs).

[↑ Back to top](#table-of-content)

---

## 11. Contig-redundancy screen

Scripts: [`contig_redundancy.slurm.sh`](scripts/contig_redundancy.slurm.sh) ·
[`contig_redundancy.py`](scripts/contig_redundancy.py).

Long-read assemblers (Canu, Flye, …) sometimes leave **small contigs that are really a
duplicated piece of a chromosome** — an uncollapsed repeat, a haplotig-style overlap, an
extra rDNA or organelle copy. Left in the assembly, the genes on them get counted as extra
copies and inflate the pangenome. Genomes with many more contigs than chromosomes are the
ones to check.

For every genome the SLURM script indexes the FASTA and self-aligns it
(`minimap2 -xasm5 -DP <asm> <asm>` — the same call `purge_dups` uses), then
`contig_redundancy.py` combines that with the vs-reference PAF from the
[synteny step](#9-synteny-vs-the-reference) and the GFF gene counts. One row per contig:

| class | meaning |
|---|---|
| `large` | ≥ 1 Mb (chromosome-scale) |
| `small-redundant` | ≥ 90 % of the contig aligns (≥ 95 % identity) to a longer contig of the same assembly — a duplicate |
| `small-placed` | small, aligns to the reference, not a duplicate — probably real |
| `small-unplaced` | small, doesn't map to the reference, not a duplicate — accessory content, divergent subtelomere, or contamination |

This is a two-step review: **screen** every genome, then **verify** the genes on the contigs
you're considering dropping. Nothing is removed until you edit the drop list in
[section 12](#12-remove-the-unwanted-contigs).

### 11.1 Screen

```bash
sbatch 01_Genome_Selection_QC/scripts/contig_redundancy.slurm.sh
```

**Parameters** — the `TUNING` block at the top of `contig_redundancy.slurm.sh` (or the same
names as `--flags` on `contig_redundancy.py`; `--help` lists all of them):

| var | default | change it? |
|---|---|---|
| `CHROM_MIN_MB` | `1.0` (megabases, i.e. 1 Mb) | Check this per species. Contigs ≥ this length are "chromosome-scale" and never flagged; lower it for a compact protist with sub-Mb chromosomes, raise it if the smallest chromosome is several Mb. |
| `DUP_FRAC` | `90` | Rarely. % of a contig that must be duplicated in a longer contig to be flagged `small-redundant`. Lower = more aggressive. |
| `DUP_MIN_ID` | `0.95` | Rarely. Min identity of the self-alignment blocks that count as duplication. |

Three more flags exist on `contig_redundancy.py` itself (`--min-aln`, `--ref-min-id`,
`--placed-pct`) but aren't in the `TUNING` block above — their defaults (`1000`, `0.80`, `50`)
work the same regardless of species, so you don't need to touch them; see
`contig_redundancy.py --help` if you ever do.

**Output:**

- `contig_redundancy/REPORT.md` — summary + flagged contigs.
- `contig_redundancy/<label>.contigs.tsv` — every contig.

Then, per flagged contig:

- `small-redundant`, **0 genes** — drop it (pure duplicated sequence).
- `small-redundant`, **has genes** — verify first, step 11.2.
- `small-unplaced` — inspect (BLAST vs `nt`, look at the annotations); may be genuine
  accessory genome, usually keep.

**Also scan the `pct_duplicated_in_longer_contig` column of `<label>.contigs.tsv` yourself.**
The `small-redundant` label only fires at ≥ 90 %, so a contig sitting just under it — e.g.
`small-unplaced` or `small-placed` at 75–89 % duplicated — is very likely redundant too and
will not be auto-flagged. Treat those as candidates: check their genes (step 11.2 takes an
explicit contig, `--contigs <id>`) and add them to the drop list if the genes check out.

### 11.2 Verify gene-carrying `small-redundant` contigs

Scripts: [`check_redundant_genes.slurm.sh`](scripts/check_redundant_genes.slurm.sh) ·
[`check_redundant_genes.py`](scripts/check_redundant_genes.py). For every gene on a flagged
contig it cuts the gene's genomic span (`samtools faidx`), `blastn`s it against the whole
assembly, keeps the best hit **not** on the gene's own contig, and checks whether that hit
lands on an annotated gene in the GFF.

Any argument after the project dir is passed straight to `check_redundant_genes.py`, so you
can scope the run:

```bash
sbatch .../check_redundant_genes.slurm.sh                          # every genome with a small-redundant contig
sbatch .../check_redundant_genes.slurm.sh <label>                  # just this genome
sbatch .../check_redundant_genes.slurm.sh <label> --contigs <id>   # just this contig — use for a
                                                                   #   borderline small-unplaced / small-placed
                                                                   #   one that 11.1 didn't auto-flag
sbatch .../check_redundant_genes.slurm.sh --also-unplaced          # add every small-unplaced contig too
```

**Parameters** (`check_redundant_genes.py --help`):

| flag | default | meaning |
|---|---|---|
| `[label ...]` | all flagged genomes | restrict to these genomes |
| `--contigs <id> ...` | — | restrict to these contigs (one genome); the way to test a contig 11.1 didn't flag |
| `--also-unplaced` | off | also run `small-unplaced` contigs, not only `small-redundant` |
| `--min-id` | `98` | % identity for a gene's off-contig hit to count as a duplicate |
| `--min-cov` | `90` | % of the gene's length that hit must cover |

**Output,** per genome, in `contig_redundancy/`:

- `<label>.redundant_genes.tsv` — one row per gene: `contig · gene · description ·
  best_off_contig_hit · pct_gene_covered · pct_identity · hit_gene · hit_gene_description ·
  verdict`. `description` / `hit_gene_description` come from the GFF `product` / `description`
  attributes, so you can see *what* each gene is.
- `<label>.redundant_genes.fa` — the extracted gene spans (for a manual BLAST, TE-aid, …).
- `<label>.redundant_genes.blast.tsv` — the raw blastn hits.

The log gives a per-contig verdict and lists every unconfirmed gene with its description:

- **DROP-OK** — every gene has a ≥ 98 % identity / ≥ 90 % coverage copy on an annotated gene
  elsewhere → add the contig to the drop list.
- **REVIEW** — some genes unconfirmed → read those rows (their descriptions tell you whether
  it's a real orphan or just a mapping artefact) before deciding.

[↑ Back to top](#table-of-content)

---

## 12. Remove the unwanted contigs

The screens in sections [10](#10-organelle-screen) and [11](#11-contig-redundancy-screen) each
end with a report naming the contigs that should not go into the pangenome — most often
**organelle contigs** (section 10), occasionally a contaminant or a redundant duplicate contig
(section 11). This section removes those contigs from `genomes/` **and** `gff_originals/`
before step 02, so every genome carries the same nuclear-only content.

Scripts: [`drop_contigs_slurm.sh`](scripts/drop_contigs_slurm.sh) ·
[`drop_contigs.py`](scripts/drop_contigs.py).

### 12.1 Write the drop list

One `<genome_basename>  <contig_id>` per line in `01_Genome_Selection_QC/drop_contigs.txt`
(`#` comments and blank lines allowed). The contig ID is the FASTA header's first token — the
same string that is in GFF column 1 (section 5).

Treat this file as a **cumulative manifest** — the running record of every contig removed
from the dataset. Add lines as you find more (organelle, then redundant, then anything else);
don't delete old ones. Re-running is safe: a contig already gone from both the FASTA and the
GFF is skipped with a note, not an error (`--strict` restores the error, to catch typos).

Read the contigs off `organelle_screen/organelle_screen.tsv`:

- **drop** a contig when the organelle covers (almost) the whole of it — a `full`,
  `full_in_scaffold`, or `MULTIPLE_COPIES` contig is the organelle itself;
- **keep** a contig where the organelle is only a fragment sitting inside a much larger
  nuclear sequence (`fragment` / a NUMT / NUPT) — that is real nuclear DNA.

...and off `contig_redundancy/REPORT.md` ([section 11](#11-contig-redundancy-screen)):

- **drop** a `small-redundant` contig once you've confirmed it carries no unique genes;
- **keep** `small-placed` / `small-unplaced` unless a BLAST check says otherwise.

```
# genome            contig to drop
AF293               mito_A_fumigatus_Af293
GCA_040126065.1     CP146368.1
GCA_040142845.1     JBBFMX010000010.1
GCA_040142845.1     JBBFMX010000011.1
GCA_040142845.1     JBBFMX010000012.1
```

### 12.2 Run it

The wrapper prints a **dry run** first — every contig it would remove, with its length (and a
note on any already gone, e.g. from an earlier run) — then **applies it in the same run**, no
separate step needed:

```bash
sbatch 01_Genome_Selection_QC/scripts/drop_contigs_slurm.sh
```

To preview only, without writing anything, call the underlying script yourself instead —
`--apply` is what turns the preview into a real write, so leave it off for the preview and add
it only once you've checked the output:

```bash
python3 01_Genome_Selection_QC/scripts/drop_contigs.py 01_Genome_Selection_QC/drop_contigs.txt              # preview only
python3 01_Genome_Selection_QC/scripts/drop_contigs.py 01_Genome_Selection_QC/drop_contigs.txt --apply       # then, once satisfied
```

| argument | meaning |
|---|---|
| `drop_contigs.txt` | the `<genome> <contig>` list (positional, required) |
| `--apply` | actually rewrite the files (default: dry run, report only) |
| `--genomes-dir DIR` | assemblies as `<genome>.fna` / `.fasta` / `.fa` (default `genomes`) |
| `--gff-dir DIR` | annotations as `<genome>.gff3` / `.gff`; `""` to skip GFFs (default `gff_originals`) |

For each genome in the list it rewrites `genomes/<genome>.fna` without those sequences and
`gff_originals/<genome>.gff3` without their features (annotation lines, `##sequence-region`
pragmas, and any `##FASTA` records) — both in the same run. Genomes not named in the list are
left untouched. A contig that is already gone from the FASTA but still present in the GFF is
fine: the script just cleans the GFF, so you can also run it after stripping the assemblies
some other way.


[↑ Back to top](#table-of-content)

---

[Back to top level](../README.md)
