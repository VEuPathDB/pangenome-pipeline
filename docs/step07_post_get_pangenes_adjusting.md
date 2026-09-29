# 07 – Post-GET_PANGENES Adjusting

> **Status: work in progress.** This step has no fixed procedure yet. This README lists what still needs to be explored
> and the scripts we have so far.

***Purpose***: Look at the pangenes that `get_pangenes` returned, and decide whether some need to be adjusted before
the results are interpreted. The plots of [step 08](step08_pangenome_visualisation.md) show where the problems
seem to be; this step is where they are examined one by one.

## Table of contents

1. [Scripts so far](#scripts-so-far)
2. [What needs to be explored](#what-needs-to-be-explored)

## Scripts so far

**`explore_multigene_pangenes.py`** (wrapper: `explore_multigene_pangenes_job.sh`). Lists the pangenes that hold more than
one gene of the same genome. For each pair of genes it reports whether their CDS overlap, the distance between them, whether
they are on the same strand, the source layer of each gene, its `intersect=` tag and whether it has a protein domain. It
needs the 7b GFFs (`intersect=`) and writes one row per gene pair. Submit it from the project folder:
`sbatch --export=ALL,LABEL=Toxo7b $PIPELINE_DIR/workflow/step07_post_get_pangenes_adjusting/explore_multigene_pangenes_job.sh`.
It reads `$STEP06_OUT/for_get_pangenes_7b_pangenes/` and `$STEP05_OUT/tagged_annotation/`, and writes to
`$STEP07_OUT/multigene_pangenes/` — by default `output/step07_post_get_pangenes_adjusting/multigene_pangenes/` in the
project folder (see [outputs](../README.md#outputs)).

**`explore_pseudogenes.py`** (no SLURM wrapper yet). Flags pseudogene-like genes by two marks kept separate: the pseudogene
biotype of the source annotation, and a broken ORF (`valid_ORFs=0`) in a Liftoff transfer. It reports them per genome and
per pangene, and writes one row for each pangene with at least one flagged gene. Genes taken from a genome's own annotation
have no ORF result, so only the biotype can flag them.

## What needs to be explored

**1. Pangenes with more than one gene in a genome.** `get_pangenes` often joins two adjacent genes of the same genome
into one pangene. In *T. gondii* they are 1,501 pangenes without filtering (7b) and 1,374 with filtering (7a). They may be
tandem duplicates, one gene annotated in pieces, or different genes that were clustered together. To find out:
- whether the genes are neighbours, overlap, or are far apart on the chromosome;
- whether they carry the same domain or complementary ones;
- whether some should be counted as separate pangenes.

**2. Invalid and pseudogene-like genes.** In *T. gondii* (7b) about two thirds of the cloud pangenes have mostly invalid
genes. It needs to be decided whether such pangenes are flagged, kept apart or removed. Genes from a genome's own annotation
have no ORF result; checking their ORFs would make the comparison fair.

**3. Genes that overlap another gene without evidence.** Many core pangenes fall in this group in *T. gondii*. One
possibility is that they are conserved only because of the real gene they overlap. It could be tested with the overlap
partner of each gene (the `<genome>_intersect_pairs.tsv` files of step 05) and how much of the gene the overlap covers.

**4. Genes annotated in pieces.** A very long gene split into several short ones would appear as many cloud pangenes. The
dynein heavy chain in *T. gondii* may be such a case, but the positions on the chromosome have not been checked. A check
for fragments (neighbouring genes of the same family, on the same strand, with no gap) could be added here.

**5. Organelle genes and nuclear insertions.** Cytochrome b domains are over-represented in the *T. gondii* cloud. Some of
them sit on a contig that looks like a repeated mitochondrial genome. Aligning these genes to the organelle genomes would show
which are organelle sequence, which are insertions in the nuclear genome, and which are poor transfers.

**6. Redundant contigs.** Genes on contigs duplicated in the same assembly count as extra copies. For *T. gondii* the
genes of the redundant contigs were removed after `get_pangenes` (step 09 of the species report). The same check is still
to be done for *A. fumigatus*, and `get_pangenes` could be run again to make the neighbourhoods exact.

**7. Transposable elements.** Genes that come from TEs and were not removed in [step 02](step02_te_annotation_clean_up.md)
should be handled here.

**8. One table per pangene.** A table with the properties of every pangene (occupancy, genes per genome, domains, invalid
share, overlap, contig size) would let the pangenes be filtered in one place instead of by separate scripts.

[↑ Back to top](#07--post-get_pangenes-adjusting) &nbsp;|&nbsp; [Back to top level](../README.md)
