# 08 – Final GFF

**Not built yet.** This page states the goal of the step, what it has to do, and the decisions still to be made.

## Goal

Steps 01–07 shape the annotation for one purpose: comparing protein-coding genes across strains. To do that, they
keep one transcript per gene, drop UTRs and non-coding genes, leave organelle genes and redundant contigs out, and
carry working tags (`basket=`, `model_from=`, `orf_status=`, `pangenome_excluded=`...).

Step 08 makes the final annotation: one conventional GFF3 per strain. It keeps the gene models chosen by the
pipeline and puts back what was left out for the pangenome. The protein domains are recomputed on the final proteins.

## Input

- `pangenome_final/gff/` (step 07 part 4): the step 05 models minus the genes removed in curation; the genes left out
  of the pangenome (e.g. mito) are still there, tagged `pangenome_excluded=<reason>`.
- `pangenome_final/pangene_ids.tsv` and the final matrix: the pangene of every gene (`<prefix>_pan_00001`...).
- Each strain's original annotation (`<step 05>/input/native/`), for what has to be put back.
- The genome FASTA of each strain.

## What the step has to do

1. **Choose the model of each gene.** Every gene has the model chosen by step 05, and some of these have a broken ORF.
   For each such gene, keep it, replace it, or drop it (see the first open decision).
2. **Put back what the pangenome left out** of the original annotation:
   - UTRs;
   - the other transcripts of each gene (step 05 keeps only the one with the longest CDS);
   - non-protein-coding genes (ncRNA, tRNA, rRNA...);
   - organelle genes (e.g. *T. gondii* apicoplast);
   - genes on the redundant contigs removed in step 01 (to discuss).
3. **Handle the curation results of step 07:**
   - genes removed in curation (`removed_genes.tsv`): out;
   - genes left out of the pangenome only (`excluded_from_pangenome.tsv`, e.g. mito): in or out;
   - genes not kept in the reference annotation (`reference_excluded.tsv`): out of the reference GFF only.
4. **Write a conventional GFF3:** gene / mRNA / exon / CDS / UTR features, valid against a GFF3 validator, with every
   attribute of the source annotation kept (Dbxref, description...). The working tags of steps 05–07 are dropped, or
   reduced to the few that a user needs.
5. **Add the pangene of each gene** as an attribute (`pangene=<prefix>_pan_00001`), so the GFF and the pangenome
   matrix can be joined.
6. **Run InterProScan on the final proteins** of all genes and write the domains from that run. Until then, a gene's
   domains are those of the protein its model came from (`model_from`), not of its own protein.
7. **Write the sequences:** protein and CDS FASTA of the final models.
8. **Check the result:** per strain, gene counts and how each gene relates to the original annotation (identical,
   updated, joined, added, removed), and the overlaps left.

## Decisions still to make

1. **Broken ORFs.** What to do with a gene whose final model has a broken ORF:
   - keep the original model when the original has a valid ORF (example asked: "GFF with only valid genes; where the
     model is invalid, keep the original");
   - keep it, tagged as a pseudogene or a broken ORF;
   - drop it;
   - repair it manually first (step 05 already repaired what miniprot could).
2. **How many GFFs.** One GFF, or the three discussed on 2026-11-10:
   - **official**: the original annotation, updated only where the pipeline's model is clearly better (a broken model
     replaced by a complete one, pieces joined into one gene);
   - **everything**: all the pipeline's models, coloured by what happened to them;
   - **pangenome**: the models used in the pangenome.
3. **Gene IDs.** Keep the original IDs (`TGME49_...`, `AFUA_...`) where there is one, or give new IDs to every gene.
   For updated or joined models, how to name them and how to keep the original ID(s) (e.g. an attribute listing the
   genes a model came from). New genes added by the pipeline need IDs either way.
4. **Genes left out of the pangenome only** (e.g. mito copies in the nuclear genome): kept in the final GFF, tagged,
   or removed.
5. **New genes added by the pipeline with a complete ORF:** add them to the official GFF or not.
6. **Genes with more than one transcript** (e.g. ME49 TGME49_211420, SR3): step 05 kept only the transcript with the
   longest CDS. Put the other transcripts back or not, and whether their CDS must also be checked (valid ORF).
7. **Genes on redundant contigs:** put back or not (they were left out because they would inflate the cloud and are probably assembly artefacts).
8. **Genes in no pangene, left out by `get_pangenes`** (CDS under about 100 bp; *T. gondii* 18, *A. fumigatus* 68,
   e.g. `TGME49_500178` in all 15 strains): kept in the final GFF or not, and with which tag. Some also have an invalid
   ORF (*T. gondii* COUG `Tg4g00431`, CAST `Tg3g02620`, p89 `Tg9g04979`) and are in no group of the broken-ORF review.
9. **Which working tags stay** in the final GFF (e.g. `orf_status`, `pangene`, the model source).
10. **A change table** for the database: one line per gene that differs from the original annotation, and why.
