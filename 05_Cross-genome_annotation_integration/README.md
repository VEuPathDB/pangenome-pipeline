# 05 – Cross-genome annotation integration (Annotation LiftOff all-vs-all)

***Purpose:*** To obtain accurate and consistent gene annotations across all genome assemblies, reducing annotation-related biases in the pangenome analysis.


A key premise of pangenome analysis is that gene annotation has to be **consistent across
every genome** in the set. If a gene is genuinely present in a genome, its annotation should
call it - a real locus missing from the pangenome matrix just because one genome's annotation
happened to miss it is an artefact, not biology. In practice, genome annotations
are frequently produced by different research groups using different tools and methodologies,
and depending on the tool/evidence used, the same real gene can get called in one strain's
annotation and missed in another's - not because the gene is truly absent, but because of how
it was annotated. 

We therefore developed a pipeline to **transfer and integrate gene
annotations**, producing a more consistent annotation across genomes of the same species.

Nevertheless, this approach has important limitations. Even when a transferred gene model contains a valid open reading frame, the pipeline cannot determine whether the gene is active or functional. A gene may have its expression affected by mutations in regulatory regions that are not evaluated by the pipeline. Moreover, pseudogenes are also considered in the pipeline; although disrupted genes may contain a genuine loss-of-function mutation, the apparently disruption could also be result from a sequencing or assembly error. Distinguishing among these possibilities is not trivial and is beyond the scope of this pipeline at this moment. Therefore, the presence of an annotated gene model should not necessarily be interpreted as evidence that the gene is functional.

This pipeline builds that consistent annotation around an **all-vs-all Liftoff transfer**:
every genome's genes get transferred onto every other genome. The **reference** - the one
genome whose annotation you trust the most, picked before anything else runs - is transferred
first as the trusted `main` layer. Then **two distinct
approaches exist to build the final annotation**, complementing each
genome's annotation with its own original annotation or the ones coming from other strains' genomes. In summary, the approaches are:

- [Step 5](#step-5-build-each-genomes-final-annotation-priority-order) - a fixed priority order
  of sources: only the reference's "good" transfers first (valid gene model, high enough
  coverage and sequence identity), then that genome's own original annotation, then every other
  genome's transfers - "good" ones before low-confidence leftovers - so one genome's weak
  transfer never blocks a better version of the same gene coming from a different genome.
- [Step 8](#step-8-build-each-genomes-final-annotation-best-model) - SADIK TO WRITE HERE




## Table of contents

1. [Native annotation check](#step-1-native-annotation-check)
2. [Remove specific genes from the GFFs (optional, manual kill list)](#step-2-remove-specific-genes-from-the-gffs-optional-manual-kill-list)
3. [Prepare the GFFs for Liftoff](#step-3-prepare-the-gffs-for-liftoff)
4. [Liftoff](#step-4-liftoff)
5. [Build each genome's final annotation (priority order)](#step-5-build-each-genomes-final-annotation-priority-order)
6. [Built annotation check (if using Step 5)](#step-6-built-annotation-check-if-using-step-5)
7. [Final annotation (if using Step 5, pre-get_pangenes)](#step-7-final-annotation-if-using-step-5-pre-get_pangenes)
   - [7a. Filter (drop unsupported genes before get_pangenes)](#7a-filter-drop-unsupported-genes-before-get_pangenes)
   - [7b. No filtering (tag, decide after get_pangenes)](#7b-no-filtering-tag-decide-after-get_pangenes)
8. [Build each genome's final annotation (best model)](#step-8-build-each-genomes-final-annotation-best-model)


## Steps


### Step 1: Native annotation check

**In short:**
- Reads every genome's own GFF, **before Liftoff runs** - each genome's annotation exactly as
  its own source produced it (plus InterPro/PFAM evidence added by step 04).
- Classifies every gene into protein-coding / pseudogene / tRNA / rRNA / other non-coding /
  unclassified, whatever naming convention the source annotation used.
- Checks how many protein-coding genes carry domain evidence.
- Counts overlapping gene pairs per genome, same-strand vs. opposite-strand, at three levels of
  strictness, plus a minor/substantial/contained severity breakdown.
- Read-only - never edits any GFF.

Scripts: [`scripts/analyze_input_gffs_job.sh`](scripts/analyze_input_gffs_job.sh) ·
[`scripts/analyze_input_gffs.py`](scripts/analyze_input_gffs.py). Read-only - it never edits any GFF.
Run once, after [step 04 - InterProScan](../04_InterProScan/README.md) - the domain-evidence
columns and plots (`has_domain`) only mean anything once InterPro/PFAM hits are written into the
GFFs.

**How genes are classified.** Different sources describe genes differently - one might call a
protein-coding gene `protein_coding_gene` directly, another plain `gene` with a
`gene_biotype=protein_coding` attribute, another neither - so this also doubles as a check that
this project's downstream scripts are actually reading each genome's gene calls correctly. For
each gene-level feature (any feature type containing "gene" - `gene`, `protein_coding_gene`,
`pseudogene`, `tRNA_gene`, ...), in order:
1. **By feature type**, if it's one of the specific ones above (`protein_coding_gene`,
   `pseudogene`, `tRNA_gene`, `rRNA_gene`, `ncRNA_gene`, ...) - no ambiguity, done.
2. Otherwise (plain `gene`, or anything else not in that list), **by attribute** -
   `gene_biotype=` (NCBI/GenBank-style) or `biotype=` (Ensembl-style).
3. Otherwise, **by its own children** - if it has a `CDS` descendant anywhere below it, it's
   protein-coding; a `tRNA`/`rRNA` descendant makes it that; a `pseudogenic_transcript`
   descendant makes it a pseudogene; other non-coding RNA descendant types
   (`ncRNA`/`lnc_RNA`/`snRNA`/`snoRNA`/`miRNA`/`misc_RNA`) bucket together as "other
   non-coding."
4. If none of the above gives an answer, the gene is **unclassified** - flagged so it gets a
   manual look rather than being silently miscounted as something it might not be.

Every distinct (genome, feature type, how it was classified) combination found is written to
`feature_attribute_table.tsv` - check this first if a genome's numbers look off, or if any
gene came back unclassified. A genome reporting zero gene-level features at all means its GFF
doesn't use "gene" anywhere in its feature type names - a convention this script doesn't
recognize, worth a manual look at that genome's GFF directly.

**Domain evidence** is the same check the filter ([step 7](#step-7-final-annotation-if-using-step-5-pre-get_pangenes)) itself uses: a protein-coding
gene "has a domain" if any of its CDS lines carries an `InterPro:`/`PFAM:` hit in `Dbxref=`.
Because this check reads step 04's output, that evidence is this project's own InterProScan
result for every genome. Genes without a hit are not necessarily wrong - many real genes have no
recognisable domain - but the numbers give a baseline for the "protein domain" evidence that
the filter step relies on, and `has_domain` in `cds_overlapping_pairs.tsv` helps decide which
gene of an overlapping pair is the more credible one.

**Overlapping genes** - protein-coding and pseudogene genes only (overlap between two
non-coding RNA genes, or between a non-coding and a coding one, isn't the
same-locus-consistency problem this pipeline is checking for) - counted per genome, split by
whether the pair is on the same strand or opposite strands, at three levels:

- **gene overlap** - do the two genes overlap at all, using the raw `gene` line's own
  start/end (this includes any annotated UTR);
- **core-span overlap** - same question, but using only the extent of each gene's CDS/exon
  (UTR excluded - ignores introns within either gene, just their outer boundary);
- **CDS-to-CDS collision** - do their actual CDS segments overlap.


**CDS-overlap severity.** For every CDS-collision pair, severity is based on how much of the
overlap covers the ***smaller*** gene's own CDS:

- **minor** (< 20%): likely just gene-model boundary noise, not necessarily a wrong call;
- **substantial** (20-80%): ambiguous, worth a look;
- **contained** (≥ 80%): one gene is essentially sitting inside the other - a strong candidate
  for one of the two being a spurious/duplicate call.

Same rule and thresholds as step 05's own `plot_final_annotation_qc.py` (post-Liftoff build), so
the two are directly comparable.

#### Running it

**Edit the settings at the top of `analyze_input_gffs_job.sh`**:

| Variable | Default | Meaning |
|---|---|---|
| `GENOMES` | `genome_order.txt` | Comma-separated genome names, or a path to a file with one per line (e.g. `genome_order.txt`) |
| `GFF_DIR` | `gff_with_interpro` | Directory with every genome's `<genome>.gff3`/`.gff` including InterPro evidence |
| `OUT_DIR` | `annotation_check` | Where to write the tables and plots |

**Run:**
```bash
sbatch analyze_input_gffs_job.sh
```

**Output**, in `OUT_DIR`:

| File | Shows |
|---|---|
| `feature_attribute_table.tsv` | every distinct (genome, feature type, how it was classified) combination found, with a count - start here if a genome's numbers look wrong |
| `gene_biotype_counts.tsv` | per genome: protein-coding (split by domain evidence present/absent), pseudogene, tRNA, rRNA, other non-coding, unclassified |
| `overlap_counts.tsv` | per genome: overlapping protein-coding/pseudogene pairs, split same-strand/opposite-strand, for the raw gene-overlap count, the core-span overlap count, and the CDS-collision count - plus, for the CDS-collision pairs, a minor/substantial/contained severity breakdown (both strands combined) |
| `<genome>_gene_classification.tsv` | one row per gene - id, feature type, category, how it was classified, domain evidence, coordinates |
| `cds_overlapping_pairs.tsv` | one row per CDS-overlapping gene pair, across every genome - both genes' id, description, coordinates, CDS length and domain evidence (has_domain), plus the overlap in bp and its severity |
| `01_gene_counts_protein_coding.png` | protein-coding genes (split: has domain evidence / doesn't) + pseudogenes, per genome, with totals labeled on top of each bar |
| `02_gene_biotype_breakdown.png` | every gene category found, per genome, totals on top |
| `03_gene_overlap_counts.png` | any overlap at all (raw `gene` line span, UTR included), per genome, split same-strand/opposite-strand, totals on top |
| `04_core_overlap_counts.png` | core-span (CDS/exon extent, UTR excluded) overlapping gene pairs, per genome, split same-strand/opposite-strand, totals on top - a subset of 03's counts, usually much smaller |
| `05_cds_overlap_counts.png` | actual CDS-to-CDS collisions, per genome, split same-strand/opposite-strand, totals on top - a subset of 04's counts, shown as its own plot since stacking all three together would make the smaller ones invisible |
| `06_cds_overlap_severity.png` | the same CDS-collision pairs as 05, both strands combined, split minor/substantial/contained - the pre-Liftoff baseline to compare against step 05's later post-Liftoff `plot_final_annotation_qc.py` output |

`cds_overlapping_pairs.tsv`'s `gene*_description` comes from `description=` on the gene line
itself (FungiDB-style) or, if that's absent, `product=` on a descendant CDS/mRNA (NCBI-style) -
whichever this genome's source annotation actually used. `gene*_has_domain` is blank for a
pseudogene (not tracked for that category, same as `gene_biotype_counts.tsv`) - for a
protein-coding gene it can be `False` even when the description itself reads like a real
function (e.g. "Has domain(s) with predicted role in..."): that description is this genome's
own prior functional annotation, while `has_domain` here is strictly whether this project's own
InterProScan run (or the original annotation) left an `InterPro:`/`PFAM:` hit in this gene's
`Dbxref=` - see the "Known limitation" further down for why a genome's original annotation
alone can show 0% (which is why this check runs on step 04's output).

### Step 2: Remove specific genes from the GFFs (optional, manual kill list)

If some genes are confirmed to be wrong and must be deleted (e.g. 5 genes that are known to be wrong in *T. gondii* M49 last annotation), strip them out before Liftoff.

Scripts: [`scripts/remove_genes_from_gff_job.sh`](scripts/remove_genes_from_gff_job.sh) ·
[`scripts/remove_genes_from_gff.py`](scripts/remove_genes_from_gff.py)

- Removing a gene removes every feature nested under it too (mRNA, exon, CDS, UTR, ...), so
  no orphaned child lines are left behind.
- The input GFFs are never edited. Every genome in the genome list gets a copy in `OUT_DIR`,
  with the kill-listed genes stripped from it (a genome with nothing on the list is copied
  unchanged), so `OUT_DIR` is a complete drop-in replacement for the input directory - point
  the following steps at it.

**Kill list**: a TSV with no header, one gene per line, two columns - `genome<TAB>gene_id`
(`#` lines and blank lines are ignored). The ID has to match the gene's `ID=` attribute in
that genome's GFF; a gene ID only applies to the genome named on its own line.

```
# kill_list_manual.tsv
ME49	TGME49_242250
ME49	TGME49_224780
```

**Edit the settings at the top of `remove_genes_from_gff_job.sh`**:

| Variable | Default | Meaning |
|---|---|---|
| `GENOMES` | `genome_order.txt` | Comma-separated genome names, or a path to a file with one per line (e.g. `genome_order.txt`) |
| `GFF_DIR` | `gff_with_interpro` | Directory with every genome's GFF to remove genes from |
| `KILL_LIST` | `kill_list_manual.tsv` | The kill list TSV described above |
| `OUT_DIR` | `gff_curated` | Where to write the curated GFFs - must differ from `GFF_DIR` |

**Run:**
```bash
sbatch 05_Cross-genome_annotation_integration/scripts/remove_genes_from_gff_job.sh
```

**Output** - everything goes to `OUT_DIR` (default **`gff_curated/`**, next to `gff_with_interpro/`, in the
project folder):

| File | What it is |
|---|---|
| `gff_curated/<genome>.gff3` | one GFF **per genome in the list** - the kill-listed genes (and everything nested under them) removed; genomes with nothing on the kill list are identical copies. Same file name as the input, so it is a drop-in replacement for `gff_with_interpro/` |
| `gff_curated/removed_genes_report.tsv` | one row per removed gene: genome, gene ID, feature type, number of GFF lines removed |
| `remove_genes_from_gff_<jobid>.out` / `.err` | the job logs, in the folder you submitted from (not in `gff_curated/`) |

`gff_with_interpro/` itself is never changed. [Step 3](#step-3-prepare-the-gffs-for-liftoff) reads
`gff_curated/`; if you skip this step (no genes to remove), point step 3 at `gff_with_interpro/` instead.

**Checking it worked:** the job's `.err` log has one block per genome with kill-listed genes -
how many genes and GFF lines were removed, per gene - and a `WARNING` for any ID that wasn't
found (a typo in the ID, or the wrong genome name on that line; nothing is removed for it).
`<OUT_DIR>/removed_genes_report.tsv` has the same numbers as a table.

### Step 3: Prepare the GFFs for Liftoff

Scripts: [`scripts/prepare_for_liftoff_job.sh`](scripts/prepare_for_liftoff_job.sh) ·
[`scripts/prepare_for_liftoff.py`](scripts/prepare_for_liftoff.py)

Two things, done on a copy - the input GFFs are never edited.

**1. Choose what Liftoff transfers.** Liftoff only transfers top-level features whose type is in
its list: `gene` by default, plus whatever a `-f` types file adds - the match is exact, so
`protein_coding_gene`, `pseudogene` or `ncRNA_gene` are not transferred unless listed.

Different sources name genes differently (FungiDB: `protein_coding_gene`, `ncRNA_gene`; NCBI:
plain `gene` for everything, told apart only by `gene_biotype=`). So each gene is first
classified - same rules as the [native annotation check](#step-1-native-annotation-check): by
feature type, else `gene_biotype=`/`biotype=`/`ebi_biotype=`, else by what's nested under it -
then retyped:

- **Protein-coding genes and pseudogenes** get their feature type (column 3) rewritten to
  `protein_coding_gene`, which is then written into `liftoff_feature_types.txt` - the list
  passed to Liftoff with `-f`. That's what makes these genes the ones Liftoff transfers.
- **Non-coding RNA genes** → `tRNA_gene` / `rRNA_gene` / `ncRNA_gene`. Not in
  `liftoff_feature_types.txt`, so Liftoff skips them - a genome's tRNAs and other non-coding
  genes are **not copied onto other genomes**. They aren't lost: they stay in that genome's own
  prepared GFF, and [step 5](#step-5-build-each-genomes-final-annotation-priority-order) puts
  them back into that genome's final annotation.

This retyping only affects genes typed plain `gene` (NCBI-style) - a gene already typed
`ncRNA_gene` (FungiDB-style) is left as it is, and a gene that can't be classified at all is
also left as it is and reported (if it's typed `gene`, Liftoff's default would still transfer it).

Example, an NCBI-style tRNA:

| | Feature type (column 3) | Attributes |
|---|---|---|
| before | `gene` | `ID=gene-X;gene_biotype=tRNA;...` |
| after | `tRNA_gene` | unchanged |

Only the feature type changes - every attribute, including the biotype, is kept.

**2. Remove the UTRs.** To simplify the GFF before Liftoff, UTR regions are excluded - every
gene is reduced to just its coding sequence, taken directly from each transcript's `CDS` lines.
This works whether an annotation lists UTRs explicitly or leaves them implicit (an exon reaching
past the CDS). For every transcript that has CDS lines:

- UTR lines under it are deleted;
- the transcript (mRNA) is shrunk to the CDS extent;
- each exon is clipped to the CDS extent, and an exon entirely outside it (a UTR-only exon) is
  deleted; exons inside the CDS extent, and all CDS lines, are untouched;
- each gene is then shrunk to cover its transcripts.

Genes without a CDS (non-coding RNA, pseudogenes with no CDS) are not touched by this part.

**Edit the settings at the top of `prepare_for_liftoff_job.sh`**:

| Variable | Default | Meaning |
|---|---|---|
| `GENOMES` | `genome_order.txt` | Comma-separated genome names, or a path to a file with one per line |
| `GFF_DIR` | `gff_curated` | GFFs to prepare - the kill-list output (step 2), or `gff_with_interpro` if step 2 was skipped |
| `OUT_DIR` | `gff_for_liftoff` | Where to write the prepared GFFs - must differ from `GFF_DIR` |

**Run:**
```bash
sbatch 05_Cross-genome_annotation_integration/scripts/prepare_for_liftoff_job.sh
```

**Output**, in `OUT_DIR`:

- one prepared GFF per genome;
- `liftoff_feature_types.txt` - the type to transfer;
- `prepare_for_liftoff_report.tsv` - per genome: UTR lines/exons/genes affected, and how many
  genes of each type were retyped.

### Step 4: Liftoff

Lifts every genome's own annotation onto every other genome's FASTA, all-against-all.

Script: [`scripts/liftoff_all_vs_all.sh`](scripts/liftoff_all_vs_all.sh). It runs
[Liftoff](https://github.com/agshumate/Liftoff) once per source-target pair, skipping any pair
whose result already exists, so it is safe to resubmit.

**Edit the settings at the top of `liftoff_all_vs_all.sh`**:

| Variable | Default | Meaning |
|---|---|---|
| `GENOMES` | `genome_order.txt` | Comma-separated genome names, or a path to a file with one per line |
| `GFF_DIR` | `gff_for_liftoff` | The prepared GFFs - [step 3](#step-3-prepare-the-gffs-for-liftoff)'s output |
| `FEATURE_TYPES` | `gff_for_liftoff/liftoff_feature_types.txt` | Which types Liftoff transfers - also step 3's output |
| `FASTA_DIR` | `genomes` | Every genome's own FASTA |
| `OUT_DIR` | `liftoff_output` | Where to write the pair results |

**LiftOff parameters**

The Liftoff call itself also has these parameters, hardcoded in the script rather than exposed
as settings above - edit them directly in `liftoff_all_vs_all.sh` if you want to change one,
otherwise every run uses:

| Flag | Liftoff's own default | Our value | What it does |
|---|---|---|---|
| `-p` | `1` | `16` | Threads to use |
| `-polish` | not included | always included | **Do not remove this one** - it is essential to get more genes with valid ORFs |
| `-copies` | not included | always included | Also search the rest of the target genome for additional copies of each gene - this is what produces `extra_copy_number`; can find real new genes, but adds extra ones that need checking |
| `-sc` (minimum identity for `-copies` to count a copy) | `1.0` | `0.9` | How similar an extra copy must be to the original to be reported |
| `-flank` | `0.0` | `0.1` | Amount of flanking sequence to align, as a fraction of gene length. It can improve gene alignment where gene structure differs between target and reference |
| `-mismatch` | `2` | `1` | Allows more variation in size and sequence - useful for divergent genomes, but not essential for most |
| `-gap_open` | `2` | `1` | same as above |
| `-gap_extend` | `1` | `1` | unchanged from Liftoff's own default |
| `-f liftoff_feature_types.txt` | `gene` only | `protein_coding_gene` | Which top-level feature types to transfer - written by [step 3](#step-3-prepare-the-gffs-for-liftoff), so non-coding RNA genes (given a different type there) are not transferred. **If step 3 wasn't run, you must build this file yourself** and decide which feature types to list - e.g. a previous pipeline used `gene`, `pseudogene` and `protein_coding_gene` together |
| `-dir` | one shared `intermediate_files/` | one folder per pair | Liftoff writes temporary files (some with a fixed name, e.g. `polish.sam`) here - two Liftoff runs sharing a folder would overwrite each other, so every pair gets its own (`liftoff_output/tmp/<source>_onto_<target>`), removed when the pair finishes; a folder left there means that pair failed |
___

Liftoff's own coverage/identity scores (`coverage=`, `sequence_ID=`) on the output are what
the "build final annotation" step's quality filter (`get_good_transfer_ids.py`) and the
filter step's overlap resolution both read - nothing recomputes them independently later.

**Run** from the project folder. Job logs go to `liftoff_output/logs/` (one per array task), which has
to exist before submitting - SLURM does not create it, and without it the job leaves no log. Either as one
job (sequential):
```bash
mkdir -p liftoff_output/logs
sbatch 05_Cross-genome_annotation_integration/scripts/liftoff_all_vs_all.sh
```
or one job per source genome, in parallel (much faster; N = number of genomes):
```bash
mkdir -p liftoff_output/logs
N=$(grep -c . genome_order.txt)
sbatch --array=1-${N} 05_Cross-genome_annotation_integration/scripts/liftoff_all_vs_all.sh
```

**Output**, in `liftoff_output/`:

- `polished/<source>_onto_<target>_liftoff.gff3_polished` - one per pair, the file step 5 reads
  with `PAIR_SUFFIX="_liftoff.gff3_polished"`;
- `<source>_onto_<target>_unmapped.txt` - one per pair.

### Step 5: Build each genome's final annotation (priority order)

This step, and [Step 6](#step-6-built-annotation-check-if-using-step-5) and
[Step 7](#step-7-final-annotation-if-using-step-5-pre-get_pangenes) that follow it, decide, gene locus by gene
locus, whose version of that gene to trust - by a strict priority order of sources. **Want a
different methodology - deciding by whichever transferred gene model is individually
best-supported, regardless of source? Skip straight to
[Step 8](#step-8-build-each-genomes-final-annotation-best-model)** - not written yet, so
following it end to end isn't possible yet.

Scripts: [`scripts/build_final_annotation_basic.sh`](scripts/build_final_annotation_basic.sh) ·
[`scripts/build_final_annotation.sh`](scripts/build_final_annotation.sh) ·
[`scripts/get_good_transfer_ids.py`](scripts/get_good_transfer_ids.py) ·
[`scripts/check_overlapping_genes.py`](scripts/check_overlapping_genes.py) ·
[`scripts/check_final_annotations_job.sh`](scripts/check_final_annotations_job.sh)

**What this step does.** Building one genome's final annotation means deciding, gene locus
by gene locus, whose version of that gene (if anyone's) to trust. This script does that by
adding genes in a strict priority order, one source at a time - once a locus is filled by an
earlier, more-trusted source, nothing added later can take it over or edit it. 

The **reference genome always goes first**, since by definition it's the annotation you trust
most for the whole project. For
every other genome, only the reference's genes that Liftoff could place *confidently* (we referred as **"good gene"**: a
valid ORF, ≥90% coverage, ≥90% identity, and the primary copy, not a Liftoff-detected
duplicate) become that starting layer, tagged `main`.

Next, wherever `main` left a locus blank, the genome's own **original annotation** gets a
chance to fill it in. This layer is tagged with the genome's own original source, e.g. `Genbank`.

After that, every **other** genome in your list takes a turn too, in the order you listed
them, twice over:
- first, round 1, only with their own confident transfers - so no genome's
weak evidence can occupy a locus before a stronger version from somewhere else gets a fair
shot at it - tagged `good_<that genome>`. 
- Only once *every* genome has had that first,
high-confidence turn does the process go around a second time (round 2), this time taking
each genome's *entire* transfer regardless of quality - a last resort for whatever's still
blank - tagged `rescued_<that genome>` (the reference's own leftovers go first here too).
Whatever's still blank after all of that just isn't in this genome's final annotation - no
genome, including the reference, offered anything usable for that locus.

```
for the genome being built right now ("this genome"):

  layer 1: the reference's genes
       this genome IS the reference?  -> use its own annotation, no filtering needed
       otherwise                      -> only the good Liftoff transfers
     ─────────────────────────────────  tagged "main"

  layer 2: + this genome's own original annotation, filling whatever layer 1 left blank
     (skipped for the reference - layer 1 already used its own annotation)
     ─────────────────────────────────  tagged with its own source, e.g. "Genbank"

  layer 3: + every OTHER genome, one at a time, in your chosen order -
           only its confident transfers, filling whatever's still blank
     ─────────────────────────────────  tagged "good_<other genome>"

  layer 4: + every OTHER genome again, same order - now its full, UNFILTERED
           transfer, filling whatever's still blank (reference's leftovers first)
     ─────────────────────────────────  tagged "rescued_<other genome>"

  = <this genome>_final_annotation.gff
```

The result, `<genome>_final_annotation.gff` (one per genome), has one row's `source` value
telling you exactly where every gene came from:

| `source` value | Meaning |
|---|---|
| `main` | The reference's own annotation as-is (for the reference itself), or its high-confidence transfer onto this genome |
| *(this genome's own original source, e.g. `Genbank`)* | This genome's own independent gene call - skipped for the reference, nothing to add |
| `good_<genome>` | That genome's high-confidence transfer onto this one |
| `rescued_<genome>` | That genome's full, unfiltered transfer, added only as a last resort |

Everything is done with AGAT's `complement` tool, which only *adds* genes at empty loci - it
never merges or reshapes an existing gene model.

#### Dependencies

- **AGAT**, specifically `agat_sp_filter_feature_from_keep_list.pl` and
  `agat_sp_complement_annotations.pl`. `build_final_annotation_basic.sh` activates it with
  `micromamba activate bioinfo` - change that one line if your cluster's AGAT environment has
  a different name.
- **Python 3**, standard library only - nothing to install.

#### What you need before starting

- Every genome's own GFF, in one directory: **`gff_for_liftoff/`**, the output of
  [prepare for Liftoff](#step-3-prepare-the-gffs-for-liftoff). This is the "original annotation" layer
  described above, and it is the prepared version on purpose: it still contains every gene,
  including the non-coding RNA genes (as `ncRNA_gene`/`tRNA_gene`), which Liftoff did not
  transfer, so this layer is what puts each genome's own non-coding genes back into its final
  annotation. Its gene types are already unified (`protein_coding_gene`) and the kill list is
  already applied. The trade-off: it has no UTRs, same as the Liftoff layers. (To keep UTRs
  instead, point `INPUT_DIR` at the GFFs from before the prepare step - but tRNAs typed plain
  `gene` then stay `gene`, and the get_pangenes preparation cannot tell them from real genes.)
- Your all-against-all Liftoff output (already done) - one
  `<genome>_onto_<other>_liftoff.gff3_polished` file per ordered genome pair.
- **The genome order** - decide which genome is the **reference** (the annotation you trust
  most) and what order to check the rest in. This becomes `ORDER` in
  `build_final_annotation_basic.sh` below - first name in the list = the reference. This
  order also affects the *rest*: when two lower-priority genomes could both fill the same
  gap, the one listed first wins (see above).

#### Running it

**Edit the settings at the top of `build_final_annotation_basic.sh`** (defaults are set for the
standard layout of this pipeline - `gff_for_liftoff/`, `liftoff_output/polished/`, `genome_order.txt`):

| Variable | Default | Meaning |
|---|---|---|
| `INPUT_DIR` | `gff_for_liftoff` | Directory with every genome's own GFF - the prepared GFFs (step 3) |
| `PAIRWISE_DIR` | `liftoff_output/polished` | Directory with your Liftoff output |
| `OUTPUT_DIR` | `expanded_annotation` | Where to write the final annotations |
| `ORDER` | `genome_order.txt` | Comma-separated genome names, **first one = the reference** you trust most (e.g. `"AF293,strain2,strain3"`), or a path to a file with one genome name per line |
| `ORIGINAL_SUFFIX` | `.gff3` | Filename suffix for each genome's own GFF in `INPUT_DIR` (the prepare step always writes `.gff3`) |
| `PAIR_SUFFIX` | `_liftoff.gff3_polished` | Filename suffix of the Liftoff output |
| `MIN_COVERAGE` | `0.9` | Minimum Liftoff coverage for a transfer to count as "good" (round 1) instead of falling back to a low-confidence "rescue" (round 2) |
| `MIN_IDENTITY` | `0.9` | Minimum Liftoff sequence identity for a "good" transfer, same idea |

Also check the `#SBATCH` lines at the very top (`--time` especially - building every genome
one after another can take a while with a lot of genomes) and the `micromamba activate
bioinfo` line (change the environment name if your cluster's AGAT isn't called that).

**Run** from the project folder. Job logs go to `expanded_annotation/logs/`, which has to exist before
submitting (SLURM does not create it). Either every genome in one job:
```bash
mkdir -p expanded_annotation/logs
sbatch 05_Cross-genome_annotation_integration/scripts/build_final_annotation_basic.sh
```
or one job per genome, in parallel (each genome's build is independent; N = number of genomes):
```bash
mkdir -p expanded_annotation/logs
N=$(grep -c . genome_order.txt)
sbatch --array=1-${N} 05_Cross-genome_annotation_integration/scripts/build_final_annotation_basic.sh
```

#### Checking it worked

- `<OUTPUT_DIR>/<genome>_final_annotation.gff` should exist for every genome once the job
  finishes - one independent file per genome, no merge step. To check all of them at once, run
  from the project folder:
  ```bash
  bash 05_Cross-genome_annotation_integration/scripts/check_final_annotations_job.sh
  ```
  It prints `OK`, `RUNNING` (the job is still in the queue or running) or `MISSING` per genome. For each
  missing genome it says how far the build got (main, round 1 or round 2), shows the end of that job's
  `.err` log, names the AGAT log to read, and prints the `sbatch --array=...` command that rebuilds only
  those genomes. It exits with an error if any genome is missing, and asks you to run it again if some
  are still running.
- `expanded_annotation/logs/build_final_annotation_<jobid>.out`/`.err` (from the `#SBATCH --output`/`--error`
  lines at the top; one pair per genome when run as an array) has one `=== Genome N/total: <genome> ===` block per genome as it works through
  `ORDER`, and ends with `All genomes processed.` or an `ERROR:` line.
- `<OUTPUT_DIR>/<genome>_overlaps.tsv` is a QC file listing overlapping gene pairs in that
  genome's final annotation - not an error signal, just worth a look. Ready as soon as the build
  finishes, one row per pair: strand, severity, overlap bp, and each gene's own ID, description,
  source layer, coordinates, CDS length, domain evidence (`has_domain`) and valid-ORF status
  (`valid_ORF`) - enough to judge a pair by hand without opening the GFF. `valid_ORF` is `NA`,
  not `False`, for a gene that never went through Liftoff (the reference's own `main` genes, and
  every genome's own original-annotation layer) - Liftoff never computed that flag for them, so
  `False` would wrongly read as "checked and found invalid". Same-strand pairs are listed first,
  since those "shouldn't normally happen" and are the most worth a look;
  [step 6](#step-6-built-annotation-check-if-using-step-5)'s `cds_overlapping_pairs.tsv` has the same
  information for every genome in one table.

### Step 6: Built annotation check (if using Step 5)

> Written against this step's own [Step 5](#step-5-build-each-genomes-final-annotation-priority-order) output
> specifically - its `source` tiers (`main`/`good_<genome>`/`rescued_<genome>`/own-annotation)
> drive several of the checks below. [Step 8](#step-8-build-each-genomes-final-annotation-best-model)
> (not written yet) picks each locus's source differently, so it will need its own version of
> this step once it exists, not this one.

Scripts: [`scripts/plot_final_annotation_qc_job.sh`](scripts/plot_final_annotation_qc_job.sh) ·
[`scripts/plot_final_annotation_qc.py`](scripts/plot_final_annotation_qc.py)

Once every genome has its `_final_annotation.gff`, this gives a first look at the raw build's
quality - before any filtering, which is [step 7](#step-7-final-annotation-if-using-step-5-pre-get_pangenes).
Sample output is in [`scripts/example_qc_plots/`](scripts/example_qc_plots/).

**Edit the settings at the top of `plot_final_annotation_qc_job.sh`**:

| Variable | Default | Meaning |
|---|---|---|
| `STRAINS` | `genome_order.txt` | Comma-separated genome names, or a path to a file with one per line - **first one = the reference**, same as the build step's `ORDER` |
| `FINAL_DIR` | `expanded_annotation` | The build step's `OUTPUT_DIR` |
| `GENOME_DIR` | `gff_for_liftoff` | Directory with every genome's own annotation, to compare the gene counts with |
| `GENOME_SUFFIX` | `.gff3` | Filename suffix for each genome's own GFF in `GENOME_DIR` |
| `OUT_DIR` | `final_annotation_qc_plots` | Where to write the plots and `summary.tsv` |

It needs Python 3 + matplotlib; the job activates the `plots` micromamba environment, which has both.

**Run** from the project folder. The job log goes to `final_annotation_qc_plots/logs/`, which has to
exist before submitting (SLURM does not create it):
```bash
mkdir -p final_annotation_qc_plots/logs
sbatch 05_Cross-genome_annotation_integration/scripts/plot_final_annotation_qc_job.sh
```

Produces, in `OUT_DIR` (`final_annotation_qc_plots/`):

| File | Shows |
|---|---|
| `01_gene_counts.png` | protein-coding genes: own original annotation vs. this pipeline's output, per genome |
| `02_valid_orf_rate.png` | % of genes with a valid ORF, per genome |
| `03_extra_copies.png` | extra-copy genes, split by valid ORF, per genome |
| `04_overlap_severity.png` | CDS-overlap severity (minor/substantial/contained), per genome |
| `05_source_composition.png` | where the genes came from (`main`/original/`good_*`/`rescued_*`), per genome |
| `06_overlap_strand.png` | the same pairs as `04_overlap_severity.png`, split same-strand vs opposite-strand instead of by severity |
| `07_overlap_domain_support.png` | the same pairs as `04_overlap_severity.png`, split by domain evidence on both genes / one / neither |
| `08_minor_overlap_position.png` | `minor`-severity pairs only (domain evidence not considered here), split by `flank`/`inside`/`single_exon` - see `minor_position` above |
| `summary.tsv` | the numbers behind every plot above, one genome per row, plus two columns not in any plot: `n_same_strand_overlap`/`n_opposite_strand_overlap` - the same pairs behind `04_overlap_severity.png`/`06_overlap_strand.png` - and `n_genes_in_overlap`, the number of *distinct* genes involved in at least one of those pairs (a gene in 3 pairs counts once, so this is always ≤ `2 × n_overlap_pairs`) |
| `extra_copy_genes.tsv` | every extra-copy gene behind `03_extra_copies.png` - genome, ID, description, location, source, `extra_copy_number`, `copy_num_ID` (groups every copy of the same locus), valid ORF, and Liftoff's own `coverage`/`sequence_ID` for that copy. Sorted per genome, most-duplicated locus first, so a genome standing out in the plot can be traced back to which genes are actually piling up copies |
| `cds_overlapping_pairs.tsv` | one row per CDS-overlapping gene pair behind `04_overlap_severity.png`, across every genome - same columns as [step 1](#step-1-native-annotation-check)'s own `cds_overlapping_pairs.tsv` (strand, severity, overlap bp, each gene's ID/description/coordinates/CDS length/domain evidence), plus `gene1_valid_ORF`/`gene2_valid_ORF` (`NA`, not `False`, for a gene that never went through Liftoff - see [step 5](#step-5-build-each-genomes-final-annotation-priority-order)), `gene1_source`/`gene2_source` - which layer each half of the pair actually came from (`main`, its own original annotation, `good_<genome>`, `rescued_<genome>`) - and, for `minor`-severity pairs only, `minor_position` and `mentions_cytochrome` (see below) |

**`minor_position` - where a `minor` overlap sits in the bigger gene's own exon structure, not
just how big it is.** Blank for `substantial`/`contained` pairs (only computed for `minor`, where
the position matters most - see below). One of:
- `flank` - the overlap sits on the bigger gene's own first or last CDS segment (its outermost
  exon on either side) - the overlap is near where the bigger gene's own model could legitimately
  be a bit off, or where a genuinely separate adjacent gene could sit just outside its real
  boundary. Most likely just gene-model boundary imprecision - **both genes are probably real**.
- `inside` - the overlap sits on one of the bigger gene's *internal* CDS segments, meaning the
  smaller gene's own body mostly sits inside one of the bigger gene's introns and only its CDS
  reaches into a middle exon. Much harder to explain as boundary noise - **the smaller gene is
  more likely spurious**, *except* when `mentions_cytochrome` is true (see next) - a
  cytochrome-named gene sitting inside another gene's intron looks like the same NUMT pattern
  discussed above, a plausible later insertion into an existing intron, not a bad gene call.
- `single_exon` - the bigger gene has only one CDS segment, so there is no intron to be "inside"
  of - the flank/inside distinction doesn't apply.

`mentions_cytochrome` - `true` if either gene's own `description=`/`product=` contains
"cytochrome" (case-insensitive) - flags the exception case above without having to grep the
description column by hand.
| `overlap_domain_summary.tsv` | per genome: of the pairs in `cds_overlapping_pairs.tsv`, how many have InterPro/PFAM domain evidence on both genes, on only one, or on neither (count and %), plus `n_genes_in_overlap` again for convenience |

**`04_overlap_severity.png`'s pair count is expected to be far higher than the pre-Liftoff baseline
([step 1](#step-1-native-annotation-check)'s `06_cds_overlap_severity.png`).** The severity rule
itself is identical code (`segments_overlap_bp`, the smaller gene's own CDS as the denominator,
the same 20%/80% thresholds), both strands combined the same way - so the two numbers are directly
comparable, on purpose. But what they measure is not the same *scale* of thing: the pre-Liftoff
count is one genome's own single-source annotation checked against itself, while the post-Liftoff
count is that genome's coordinate system carrying transferred copies from every other genome in the
panel, layered on top of each other. `summary.tsv`'s `n_same_strand_overlap`/`n_opposite_strand_overlap`
columns show why: AGAT's own complement step (used to build each layer in
[step 5](#step-5-build-each-genomes-final-annotation-priority-order)) only ever rejects a new gene for colliding
with a **same-strand** CDS already placed - an opposite-strand collision, or two genes added in the
same layer that collide with each other, are never blocked at that stage. Same-strand pairs are the
ones that "shouldn't normally happen" (same convention as `check_overlapping_genes.py`'s own output)
and are worth checking against `n_same_strand_overlap` per genome; a high opposite-strand count is
the expected, still-unfiltered result of combining many genomes' gene calls, and is exactly what
[step 7](#step-7-final-annotation-if-using-step-5-pre-get_pangenes)'s filter resolves next.

One thing worth knowing before reading `02_valid_orf_rate.png`: the reference genome's own
`main` genes never go through Liftoff (they're used as-is), so they never get the
`valid_ORF=` flag Liftoff computes - the reference's own bar in that plot will look very low
as a result. That's not a real quality problem, just this metric not applying to genes that
were never transferred.

**How `04_overlap_severity.png`'s minor/substantial/contained categories work.** A pair only
counts as "overlapping" here if their **CDS** (coding sequence) actually collides - two genes
whose UTRs touch, or whose introns overlap another gene's body without any CDS-to-CDS
collision, don't count. For a gene with more than one exon, every exon's CDS segment is
checked against the other gene's, and the overlap is the *total* bp shared across all of
them combined - not just whether the two genes' start-end spans overlap (a real overlap can
sit in just one exon while the rest of a multi-exon gene sits nearby without touching
anything).

The severity itself is the fraction of the **smaller** gene's own CDS length that the
overlap covers:

- **contained**: overlap ≥ 80% of the smaller gene's CDS - one gene is essentially sitting
  inside the other.
- **substantial**: 20-80%.
- **minor**: < 20%.

Sized against the smaller gene on purpose, not the bigger one - a small worked example: a
438bp overlap is only 22% of a 1986bp gene, but 100% of a 438bp one. Sizing against the
bigger gene alone would make that pair look "minor" when one of the two genes is entirely
swallowed by the overlap.

If a gene has more than one **transcript** (isoform - different from exons above), only its
longest-CDS isoform is used for this measurement, so one gene contributes one severity value
per pair, not one per isoform. In practice this barely matters for this project's genomes -
essentially all genes here have exactly one isoform.

### Step 7: Final annotation (if using Step 5, pre-get_pangenes)

> Both 7a and 7b below are written against [Step 5](#step-5-build-each-genomes-final-annotation-priority-order)'s
> output specifically - the case-based overlap rule leans on its `main`/`good_<genome>`/
> `rescued_<genome>`/own-annotation source tiers throughout. This isn't
> [Step 8](#step-8-build-each-genomes-final-annotation-best-model)'s equivalent - it has no
> source-tier concept to lean on, so its own overlap-resolution rule (not written yet) will be
> entirely different.

Two alternative approaches, both starting from [step 6](#step-6-built-annotation-check-if-using-step-5)'s
`<genome>_final_annotation.gff` files - they disagree on *when* to decide a gene is spurious.
7a decides now, before `get_pangenes` runs, using only what's knowable from one genome's own
build. 7b decides nothing yet - it tags every gene with what it overlaps and passes everything
through unfiltered, deferring the call to after `get_pangenes` has run, when a pangene's real
presence/absence pattern across the whole panel is known - much stronger evidence than anything
computable from a single genome's GFF (see 7b below for why this matters).

#### Background: what needs resolving

##### Extra copies
With the approach above, we end up transferring all genes from all genomes to every other
genome. When the parameter `-copy` is used, we might see a high number of extra copies. We need
to decide if we trust them. On this project's real *T. gondii* data, extra copies turned out to
be **genuinely enriched on short contigs by density** (genes/Mb - tens of times higher on
10-100kb contigs than on chromosome-scale ones), but that doesn't mean every extra copy on a
large contig is trustworthy either, or every one on a small contig isn't: real tandem gene-family
arrays turn up sitting on ordinary chromosome-scale contigs too (e.g. a 21-copy cluster on one
*T. gondii* strain's chromosome IX, coverage/identity ~99-100% but only 1 of the 21 copies had a
valid ORF - the rest read as pseudogenized repeat units). Contig size and domain presence on
their own don't reliably separate a real extra copy from a spurious one - a **valid ORF is the
one unambiguous signal**, so that's what this step filters on: extra-copy genes (`extra_copy_number
> 0`) are kept only if they have a valid ORF, dropped otherwise, regardless of contig size or
domain evidence.

##### Rescued-coverage floor

`rescued_<strain>` is this pipeline's lowest-confidence source tier - round 2 of
[Step 5](#step-5-build-each-genomes-final-annotation-priority-order)'s build, a strain's *entire* transfer
added unfiltered as a last resort, unlike `good_<strain>` which already had to clear a
coverage/identity bar to get in. Before anything else, every `rescued_<strain>` gene below
`--min-rescued-coverage` (default 0.5) Liftoff coverage is dropped outright - the one basic
quality floor that tier gets, independent of ORF/domain/overlap status.

##### Overlapping genes

We observed that after obtaining the final annotation, several genes are found to be
overlapping. Overlaps under a minimum bp floor (default 20bp) are treated as boundary-rounding
noise, not real overlaps, and ignored entirely - most of these turned out, on inspection of
this project's own real data, to be 1-2bp touches from Liftoff's coordinate polishing. Of what's
left, severity is scored against whichever gene's CDS is smaller (sizing against the bigger
gene alone can make a real problem look tiny): **minor** (<20% of the smaller gene's own CDS) vs.
**substantial**/**contained** (≥20%).

A **minor** overlap only gets kept automatically if it's also **flank** - it only touches the
bigger gene's outermost CDS segment, not an internal one, so the two gene models most likely
just disagree at the boundary and both genes are probably real (same `flank`/`inside`
classification as [step 6](#step-6-built-annotation-check-if-using-step-5)'s `minor_position` column). A
minor overlap that's **inside** (touches an internal CDS segment - i.e. the smaller gene mostly
sits in what's actually an intron) or **single_exon** doesn't get that pass - it's still more
likely a spurious call, so it's resolved the same as substantial/contained, see the rule below.

**Known limitation - "protein domain" is only as good as the source annotation.** Domain
evidence is detected purely by scanning a gene's own CDS lines for an InterPro/PFAM hit in
`Dbxref=` - nothing more, no external database lookup. If a genome's *original* annotation
never recorded that field at all, every gene sourced from that genome (`main`,
`good_<that genome>`, `rescued_<that genome>`) will test as having no domain,
unconditionally - not because those genes lack real domains, but because that annotation
source never recorded them. This project's own reference (AF293/FungiDB) is exactly this
case: 0 of its ~32,000 CDS lines carry an InterPro/PFAM hit. Measured effect on this
project's own 15-genome run: that genome's extra-copy genes on long contigs are always
dropped (the "valid ORF AND domain" rule in the filter step can never be satisfied), and its
genes always lose a non-`main` overlap dispute against a domain-annotated genome's gene unless
the source-tier/size tiebreak happens to save them - 16 real `rescued_AF293` genes were
dropped this way. That is why [step 04 - InterProScan](../04_InterProScan/README.md) runs
before everything else in this step: it writes real domain evidence into every genome's GFF, so
the rule works from each genome's own InterProScan result rather than whatever the source
annotation happened to record. **[Step 1's native annotation check](#step-1-native-annotation-check)**, run on
step 04's output, shows in `gene_biotype_counts.tsv` (the `protein_coding_with_domain`/
`protein_coding_without_domain` columns) how much domain evidence each genome actually has.

**Overlap-resolution rule.** Every overlapping pair falls into exactly one of six cases,
decided by whether each gene is `main` and whether it carries domain evidence:

| # | Case | Rule |
|---|---|---|
| 1 | `main` × `main` | keep both |
| 2 | `main` × other, other **has** a domain | keep both |
| 3 | `main` × other, other has **no** domain | keep only `main` |
| 4 | other × other, **exactly one** has a domain | keep only the one with a domain |
| 5 | other × other, **both** have a domain | keep both |
| 6 | other × other, **neither** has a domain | see below - no domain evidence anywhere to decide by |

Cases 1-5 need nothing beyond domain evidence (or being `main`). Case 6 - domain evidence
exhausted on both sides - falls back to how well-supported each gene is **across the whole
panel**, not just this one genome:
- a gene present in only 1 genome of the whole `--strains` panel is dropped outright - a
  single occurrence with no domain and no reference backing is the weakest possible support
  (checked against this project's real *T. gondii* data: this is also exactly the profile of
  most genuine strain-private gene calls, so this is a soft call, not a certainty);
- between two genes both present in >1 genome, keep whichever has more `good_<genome>`-sourced
  occurrences across the whole panel (a high-confidence Liftoff transfer is stronger support
  than raw presence count, which also counts low-confidence `rescued_<genome>` copies); on a
  tie, keep the one with the longer CDS.

On the *T. gondii* panel's real data (15 genomes, `--min-overlap-bp 20`): 7,777 overlapping
pairs total, of which 1,122 are case 6. Of the 914 distinct genes involved, this rule drops 672
(present in only 1 genome) and keeps 242 (present in >1, won on `good_<genome>` count or CDS
length).

**Known limitation - domain evidence on a transferred gene reflects the *original* gene, not
necessarily the transferred copy.** InterProScan (step 04) only ever runs on each genome's own
original proteins. A `good_<X>`/`rescued_<X>` gene's domain evidence is looked up from its
*source* genome's original protein, by `protein_id` - never re-checked against the transferred
copy's own sequence as it ended up in this genome. For a `good_<X>` transfer (already passed
the ≥90% coverage/identity/valid-ORF bar) that's a safe enough assumption. For a `rescued_<X>`
transfer it isn't - that's exactly the tier that failed that quality bar, so the transferred
copy's real sequence (ORF, coverage) may no longer match what the domain call was made on. For
now this is accepted as good enough for the intermediate builds; the plan is to re-run domain
detection for real, on the actual final protein set, when the final GFF is assembled for
publication (the same pass that cleans up the internal `source=` tags).

### 7a. Filter (drop unsupported genes before get_pangenes)

Scripts: [`scripts/filter_final_annotation_job.sh`](scripts/filter_final_annotation_job.sh) ·
[`scripts/filter_final_annotation.py`](scripts/filter_final_annotation.py) ·
[`scripts/plot_filter_qc_job.sh`](scripts/plot_filter_qc_job.sh) ·
[`scripts/plot_filter_qc.py`](scripts/plot_filter_qc.py)

Applies both filters above (extra copies, then overlaps - the rule table above) and writes one
filtered GFF per genome. Unlike every other script in this step, it reads **all** genomes'
`<genome>_final_annotation.gff` in one run, because case 6 above needs to know each gene's
support panel-wide, not just within its own genome - the `--strains` list has to be the full
panel, not one genome at a time. Every kept gene with no valid ORF is relabeled
`protein_coding` → `pseudogene` in its biotype attribute (not its feature type), and every kept
gene's ID is renumbered to `<id-prefix><strain-index>g<00001, ...>` in genomic order (original
ID kept in `description=`). Full rule details are in `filter_final_annotation.py`'s own header
comment.

#### Dependencies

- **Python 3**, standard library only, for filtering.
- **matplotlib** (`pip install matplotlib`) for the QC plots only.
- No AGAT, no Liftoff needed for this step - those were only for the build step above.

#### What you need before starting

- Every genome's own `<genome>_final_annotation.gff`, from the build step's `OUTPUT_DIR`.
- **The same genome order** used in the build step (`ORDER`, first = reference) - reused here
  so strain-index numbering and the reference/`main` logic stay consistent.

#### Running it

**Edit the settings at the top of `filter_final_annotation_job.sh`**:

| Variable | Default | Meaning |
|---|---|---|
| `STRAINS` | `genome_order.txt` | Comma-separated genome names or a file with one per line, **first = the reference**, same as the build step's `ORDER` |
| `FINAL_DIR` | `expanded_annotation` | Directory with every genome's raw `<genome>_final_annotation.gff` (the build step's `OUTPUT_DIR`) |
| `OUT_DIR` | `filtered_annotation` | Where to write the filtered GFFs + `_filter_report.tsv` files |
| `MIN_RESCUED_COVERAGE` | `0.5` | Every `rescued_<genome>` gene below this Liftoff coverage is dropped outright |
| `MIN_OVERLAP_BP` | `20` | CDS-to-CDS overlaps under this many bp are ignored entirely |
| `ID_PREFIX` | `Afu` | Renamed gene IDs look like `<ID_PREFIX><strain-index>g<00001, ...>` - **edit this per species**, the default only makes sense for *A. fumigatus* (matches its reference's own native gene ID style, `Afu1g00100`) |

**Run** (one single SLURM job, works through every genome in `STRAINS` in one Python process -
the cross-genome case-6 pass needs all of them loaded together anyway):
```bash
mkdir -p filtered_annotation/logs
sbatch 05_Cross-genome_annotation_integration/scripts/filter_final_annotation_job.sh
```

Then, once that's done, **edit the settings at the top of `plot_filter_qc_job.sh`**
(`STRAINS`, `ORIGINAL_DIR`, `FILTERED_DIR`, `GENOME_DIR`/`GENOME_SUFFIX`, `OUT_DIR`) and:
```bash
sbatch 05_Cross-genome_annotation_integration/scripts/plot_filter_qc_job.sh
```

**Output**, in `OUT_DIR`:

| File | Shows |
|---|---|
| `01_gene_counts.png` | total protein-coding genes, before vs. after filtering, per genome |
| `02_extra_copies.png` | extra-copy genes remaining, before vs. after, per genome |
| `03_removal_reasons.png` | why genes were removed, per genome (stacked by reason) |
| `04_overlap_severity.png` | CDS-overlap severity (minor/substantial/contained), before vs. after |
| `05_source_composition.png` | where the genes came from (`main`/original/`good_*`/`rescued_*`), before vs. after |
| `summary.tsv` | the numbers behind every plot above, one genome per row |

#### Checking it worked

- `<OUT_DIR>/<genome>_final_annotation_filtered.gff` should exist for every genome once the
  job finishes.
- `<OUT_DIR>/<genome>_filter_report.tsv` lists every removed gene and why - nothing is ever
  silently dropped.
- `<OUT_DIR>/filter_summary.tsv` - one row per genome: `n_original`/`n_kept`/`n_removed`.
- The job's `..._<jobid>.out`/`.err` (from the `#SBATCH --output`/`--error` lines at the
  top) has one `=== Genome N/total: <genome> ===` block per genome, and ends with
  `All genomes filtered. Outputs are in ...`.

Once every genome has its `_final_annotation_filtered.gff` and the QC plots look right, move
on to [step 06 - GET_PANGENES](../06_GET_PANGENES/README.md) to prepare and run
`get_pangenes` itself.

### 7b. No filtering (tag, decide after get_pangenes)

Scripts: [`scripts/tag_final_annotation_job.sh`](scripts/tag_final_annotation_job.sh) ·
[`scripts/tag_final_annotation.py`](scripts/tag_final_annotation.py)

**Why this exists.** 7a's case 6 rule (and, more subtly, cases 3/4 too) can be wrong in a
specific way: if a gene overlaps something in only one genome of the panel because that one
genome's `main`/other gene model has a locally bad Liftoff boundary - not because the gene
itself is spurious - dropping it *in that one genome* still creates a false gene-absence call
for a real, panel-wide gene, right where `get_pangenes` needs presence/absence to be accurate.
7a's panel-wide case-6 rule already guards against the worst of this (a gene present in >1
genome is never dropped just for colliding once), but the same risk exists, less mitigated, in
cases 3 and 4. 7b sidesteps the question entirely: don't decide yet, tag instead, and use
`get_pangenes`' own pangene-level presence/absence pattern - much stronger evidence than
anything computable pre-pangenome - to make the real call afterward.

Every gene in every `<genome>_final_annotation.gff` gets one new attribute:

| `intersect=` value | Meaning |
|---|---|
| `NA` | No CDS-to-CDS overlap ≥ `--min-overlap-bp` with anything |
| `with_main` | Overlaps the reference-derived `main` gene |
| `with_other_with_domain` | Overlaps a non-`main` gene that itself has a protein domain |
| `with_other_no_domain` | Overlaps a non-`main` gene with no protein domain either |

If a gene overlaps more than one other gene, the tag is the single most informative of the
above (priority `with_main` > `with_other_with_domain` > `with_other_no_domain`) - every
partner is still listed, one row per pair, in `<genome>_intersect_pairs.tsv`. Nothing is
removed, renumbered, or relabeled - every other line is byte-identical to the input.

`get_pangenes` doesn't read or preserve this attribute (it only looks at coordinates/feature
type for its WGA-based clustering), so it isn't needed to survive into `get_pangenes`' own
output - `get_pangenes`' `pangene_matrix.tab` references genes by their original gene ID, so
the tag can be joined back onto each pangene cluster afterward by ID. That join, and the
pangene-level confidence call it produces, isn't built yet - planned as part of
[step 07 - Post-GET_PANGENES adjusting](../07_Post_GET_PANGENES_adjusting/README.md).

#### Dependencies

- **Python 3**, standard library only.

#### What you need before starting

- Every genome's own `<genome>_final_annotation.gff`, from the build step's `OUTPUT_DIR`.
- **The same genome order** used in the build step (`ORDER`, first = reference).

#### Running it

**Edit the settings at the top of `tag_final_annotation_job.sh`**:

| Variable | Default | Meaning |
|---|---|---|
| `STRAINS` | `genome_order.txt` | Comma-separated genome names or a file with one per line, **first = the reference** |
| `FINAL_DIR` | `expanded_annotation` | Directory with every genome's raw `<genome>_final_annotation.gff` (the build step's `OUTPUT_DIR`) |
| `OUT_DIR` | `tagged_annotation` | Where to write the tagged GFFs + `_intersect_pairs.tsv` files |
| `MIN_OVERLAP_BP` | `20` | CDS-to-CDS overlaps under this many bp are ignored entirely |
| `ID_PREFIX` | `Afu` | Renamed gene IDs look like `<ID_PREFIX><strain-index>g<00001, ...>` - **edit this per species**, the default only makes sense for *A. fumigatus* (matches its reference's own native gene ID style, `Afu1g00100`) |

**Run:**
```bash
mkdir -p tagged_annotation/logs
sbatch 05_Cross-genome_annotation_integration/scripts/tag_final_annotation_job.sh
```

#### Checking it worked

- `<OUT_DIR>/<genome>_final_annotation_tagged.gff` should exist for every genome, same gene
  count as the input.
- `<OUT_DIR>/<genome>_intersect_pairs.tsv` lists every overlapping pair behind each tag.
- `<OUT_DIR>/tag_summary.tsv` - one row per genome: counts for `NA`/`with_main`/
  `with_other_with_domain`/`with_other_no_domain`.
- The job's log ends with `All genomes tagged (nothing removed). Outputs are in ...`.

### Step 8: Build each genome's final annotation (best model)

> ⚠️ **Not written yet** — placeholder.

[Back to top level](../README.md)
