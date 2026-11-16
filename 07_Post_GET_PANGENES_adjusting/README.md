# 07 – Post-GET_PANGENES Adjusting and Visualisation

***Purpose***: Fix what can be fixed safely in the pangenome built by `get_pangenes` ([step 06](../06_GET_PANGENES/README.md)),
look at it, remove the genes found to be false during curation, and look at it again. Only the matrix and a copy of
the GFFs change; the outputs of steps 05 and 06 are never modified.

The step has five parts, and a check of the reference before the removal:

| Part | What it does | Job | Output |
|---|---|---|---|
| 1. Inspection | tables: why genes are missing, why some pangenes hold two genes of one strain | `inspect_pangenome_job.sh` | `pangenome_inspection/` |
| 2. Safe fix | moves back together the genes that `get_pangenes` split into two pangenes; lists the pangenes that look false | `fix_pangenome_safe_job.sh` | `pangenome_fixed/` |
| 3. Preliminary look | figures and an interactive page of the fixed pangenome, to decide what to remove | `plot_all_job.sh`, `build_explorer_job.sh` with `SET=preliminary` | `pangenome_plots/preliminary/`, `pangenome_explorer/preliminary/` |
| Reference gene check | one page with the genes the pipeline added to the reference, their evidence (RNA-seq, intron junctions), and a choice per gene | `reference_check_job.sh` | `reference_check/` |
| 4. Removal after curation | removes the panel genes listed in `genes_to_remove.tsv`, and the single genes listed in `strain_genes_to_remove.tsv`, from the GFFs and the matrix; leaves the panel genes of `genes_to_exclude_from_pangenome.tsv` out of the matrix only; gives the pangenes their final IDs | `remove_panel_genes_job.sh` | `pangenome_final/` |
| 5. Final look | the same figures and page, on the final pangenome | `plot_all_job.sh`, `build_explorer_job.sh` | `pangenome_plots/final/`, `pangenome_explorer/final/` |

Parts 1 and 2 run once for each `get_pangenes` run. Parts 4 and 5 can be repeated: after each round of curation, add
the new genes to the list and run them again. Part 4 always starts from `pangenome_fixed/` and rewrites
`pangenome_final/` from scratch, so nothing is removed twice.

**User guides**, with screenshots, for the two pages people look at:
[Pangenome Explorer – user guide](explorer_guide/README.md) and
[Reference gene check – user guide](reference_check_guide/README.md). Both
pages are still under construction.

## Table of contents

- [Dependencies](#dependencies)
- [Input](#input)
- [Part 1 – Inspection](#part-1--inspection)
  - [Absences](#absences)
  - [Pangenes with more than one gene of a strain](#pangenes-with-more-than-one-gene-of-a-strain)
- [Part 2 – Safe fix](#part-2--safe-fix)
  - [Low-confidence pangenes, and genes in no pangene](#low-confidence-pangenes-and-genes-in-no-pangene)
- [Part 3 – Preliminary look](#part-3--preliminary-look)
  - [Figures](#figures)
  - [Interactive explorer](#interactive-explorer)
- [Reference gene check](#reference-gene-check)
- [Part 4 – Removal after curation](#part-4--removal-after-curation)
- [Part 5 – Final look](#part-5--final-look)
- [Running it](#running-it)
- [Outputs](#outputs)
- [Scripts](#scripts)

## Dependencies

| Program | Used by |
|---|---|
| python3 (3.9 or newer) | every script; the explorer, the inspection tables and the fixes need only the standard library |
| numpy, scipy, matplotlib | the figures (`plot_*.py`) |
| pandas, seaborn | the POCS heatmap (`plot_pocs.py`) |
| blastn, blastp, makeblastdb | the inspection (`audit_absences.py`, `inspect_multigene_pangenes.py`) |

The jobs activate the author's micromamba environment `synteny`; edit the `### ADAPT ###` block at the top of each job
to use your own.

**InterPro entry list.** The enrichment figure and the explorer name each InterPro entry from EBI's `entry.list`.
Compute nodes may not reach the internet, so download it once on the login node:
```bash
curl -L -o ~/hc-storage/interpro_entry.list https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list
```
The jobs read it from `ENTRY_LIST`; without it, the entries are labelled by their IPR accession.

## Input

| Input | From |
|---|---|
| `pangene_matrix.tr.tab`, `pangene_matrix_genes.tr.tab`, `POCS.matrix.tab` | [step 06](../06_GET_PANGENES/README.md) (`get_pangenes` writes them in a subfolder; the jobs find them by name) |
| the final GFFs (`basket_annotation/gff_prefer_panel_model_final/`) | [step 05](../05_Basket_annotation_integration/README.md): their `basket=`, `orf_status`, `model_origin` and `tandem_array_copy` tags are used |
| `winners.tsv`, `decide/absent.tsv`, `not_eligible.tsv`, `stacked_dropped.tsv` | step 05: PSAURON of every panel gene, and why step 05 left genes out |
| the genomes and the Liftoff transfers step 05 read | step 05 (inspection only) |
| `genes_to_remove.tsv`, `strain_genes_to_remove.tsv`, `genes_to_exclude_from_pangenome.tsv` | written by you during curation (part 4) |

A **panel gene** is a gene of step 05: one basket, tagged `basket=` on every gene that represents it in every strain.
The PSAURON score is per panel gene. Removals are per panel gene, or per single gene (part 4).

## Part 1 – Inspection

[`inspect_pangenome_job.sh`](scripts/inspect_pangenome_job.sh) writes two sets of tables. Nothing is changed.

### Absences

An absence is a pangene with no gene in a genome. [`audit_absences.py`](scripts/audit_absences.py) finds where the
missing gene should be (between the genes that flank it in the other genomes), looks at what that genome has there, and
gives each absence a category:

| Code | Category | Meaning |
|---|---|---|
| 1a | same_gene_movable | the same gene is there, in another pangene with no genome in common: can be moved |
| 1b | same_gene_shared | the same gene is there, in a pangene that also has genes in genomes with this pangene (paralog?) |
| 1c | same_gene_no_pangene | the same gene is there, but in no pangene (CDS < 100 bp) |
| 2a | other_gene_no_transfer | another gene is in the place, and there is no Liftoff transfer of this gene |
| 2b | other_gene_over_transfer | another gene is in the place, over a Liftoff transfer of this gene |
| 3 | transfer_not_used | nothing in the place but a Liftoff transfer of this gene, not used by step 05 |
| 4 | empty_place_shorter | nothing in the place, which is much shorter than in the other genomes: likely a deletion |
| 5 | empty_place | nothing in the place, no transfer: the gene is missing |
| 6 | contig_edge | the place is at a contig end: cannot tell |
| 7 | neighbours_split | the flanking genes are on different contigs or far apart: assembly break or rearrangement |
| 8 | region_absent | the region is not found in this genome |

"The same gene" means the DNA of the gene in place aligns with the genes of the pangene (blastn, on the CDS span and on
the spliced CDS).

### Pangenes with more than one gene of a strain

[`inspect_multigene_pangenes.py`](scripts/inspect_multigene_pangenes.py) compares the proteins of the genes that share
a cell with a reference protein of the pangene (blastp), and says why they are together: `tandem_copy` (a tandem-array
copy added by step 05), `copy`, `split` (one gene annotated in pieces), `fragment`, `unrelated` (a different gene
lumped in by `get_pangenes`), `other` or `no_reference`.

## Part 2 – Safe fix

`get_pangenes` groups genes by their place in the genome alignments. Sometimes it splits one gene into two pangenes,
with different strains in each. The gene then looks absent from both. [`fix_pangenome_safe.py`](scripts/fix_pangenome_safe.py)
moves such genes back, but only when it is safe:

- the gene sits at the expected place of the pangene and has the same sequence as its genes (absence code 1a);
- the two pangenes have no strain in common, so no strain ends up with two genes;
- the gene is claimed by one pangene only (a gene of a family that two pangenes claim is left alone);
- two pangenes that claim each other's genes merge into the one with more genomes.

Then it joins the pangenes of one gene that the moves above could not put back together, e.g. a gene that
`get_pangenes` split into three or more pangenes (each piece claims the genes of the others). Step 05's `basket=` tag says which genes are the same gene; the
pangenes holding the genes of one basket are joined into the one with the most strains, but only when:

- every one of these pangenes holds only genes of that basket;
- no strain has two genes of the basket.

In gene families a gene can carry the tag of a related panel gene, so a split where a pangene also holds genes of other
baskets, or a strain has two genes of the basket, is not joined: it is listed in `basket_splits_not_joined.tsv`, for
the curation. Every join is in `joins.tsv`.

Only the matrix changes; no gene is added to or removed from the GFFs, and no Liftoff transfer is brought back. A
pangene left empty disappears. Every absence left is explained in `absences_explained.tsv`, with a verdict:

| Verdict | Absence codes | Meaning |
|---|---|---|
| absent | 2a, 4, 5, 8 | the gene is not in that strain |
| assembly | 6, 7 | cannot tell: contig end or break |
| copy_dropped | 2b, 3 | a Liftoff copy is there, but step 05 did not use it; the column `step05` says why |
| elsewhere | 1a not moved, 1b, 1c | the gene is there, in another pangene or in none |


### Low-confidence pangenes, and genes in no pangene

After the fix, [`summarize_pangenes.py`](scripts/summarize_pangenes.py) writes two lists on the fixed pangenome, to help the curation
(nothing is removed), in `pangenome_fixed/pangene_lists/`:

- `pangene_quality.tsv`: for every pangene, the evidence to judge whether it is a real gene, sorted from the most
  invalid:
  - in how many strains its genes have an invalid ORF. A strain is valid when at least one of its genes in the pangene
    has a valid ORF, invalid when all are invalid or pseudogenes, and not judged when its ORFs are broken by the
    assembly. Pangenes with an invalid ORF in **more than 50%** of the judged strains are marked
    `low_confidence=yes` (`pseudogene` when they hold a reference pseudogene);
  - **where the ORF breaks**, read on each invalid gene's CDS: `start_lost`, `stop_lost`, `internal_stop@aaN` (first
    stop codon inside, at codon N) or `frameshift(len%3=R)`. The genes of a panel gene come from the same model, so the
    same defect in several strains (`strains_with_that_defect`, e.g. 9 of 10) is likely one shared mutation, or an error
    of the model; several different defects (`distinct_defects`) point to independent losses. A frameshift without an
    internal stop gives no position, so two different frameshifts can look the same. Needs the genomes (`GENOME_DIR`);
  - the PSAURON score of its panel gene, and its InterPro domains (`has_domain`, and each domain with the number of
    genes that carry it).

  A broken ORF with a low PSAURON score and no domain points to a false gene; a broken ORF with a high score, a domain
  and the same defect in every strain points to a real gene with a wrong model.
- `genes_not_in_pangenome.tsv`: genes of the GFFs that `get_pangenes` left out (CDS < 100 bp), with their panel gene,
  length, ORF status and PSAURON.

## Part 3 – Preliminary look

The fixed pangenome is looked at before deciding what to remove. It is looked at after the safe fix, not straight
from `get_pangenes`, so that what the fix already solves does not show up as a problem. Figures and the explorer often
show problems that no rule catches. For example, a transposable element present once in every strain, each time at a
different place, fills the cloud with single-strain pangenes and shows up as an enriched domain there. Genes found this
way are added to `genes_to_remove.tsv` (part 4). The pangenes of `pangene_quality.tsv` (part 2) can be
checked one by one in the explorer.

### Figures

[`plot_all_job.sh`](scripts/plot_all_job.sh) runs [`plot_all.sh`](scripts/plot_all.sh) and
[`plot_domain_enrichment.py`](scripts/plot_domain_enrichment.py). Every file starts with `<LABEL>_<SET>`, for example
`Toxo_preliminary_occupancy.png`. PNG only.

**Occupancy categories.** Occupancy is the number of genomes (N) a pangene is present in. The categories are cut in
numbers of genomes, not percentages, and each pangene falls in the first that fits:

| Category | Rule | Setting (default) |
|---|---|---|
| Core | missing from at most `CORE_MISSING` genomes | `CORE_MISSING` (0: present in all N) |
| Soft-core | not core, missing from at most `SOFTCORE_MISSING` genomes | `SOFTCORE_MISSING` (1: present in N-1) |
| Cloud | not core or soft-core, present in at most `CLOUD_MAX` genomes | `CLOUD_MAX` (2: present in 1 or 2) |
| Shell | everything else | |

With the defaults, 15 genomes give cloud 1-2, shell 3-13, soft-core 14, core 15; 17 genomes give cloud 1-2, shell
3-15, soft-core 16, core 17. With many genomes the defaults may be too strict. The settings are at the top of
`plot_all_job.sh` and `build_explorer_job.sh` (use the same values in both), or on the command line, for example
core = missing from at most 1 genome, soft-core = at most 3, cloud = at most 5:
```bash
sbatch --export=ALL,CORE_MISSING=1,SOFTCORE_MISSING=3,CLOUD_MAX=5 07_Post_GET_PANGENES_adjusting/scripts/plot_all_job.sh
```

| Figure | File | What it shows | What to look for |
|---|---|---|---|
| Pangenome and core size | `_pancore`, `_pancore_separate` | genomes added one at a time in 100 random orders, with fitted curves (the `.log` has the fit) | is the pangenome still growing, is the core stable |
| Pangene occupancy | `_occupancy` | one bar per occupancy, coloured by category | how much is core, shell and cloud |
| Soft-core pangenes missing per genome | `_softcore_missing` | how many soft-core pangenes each genome misses | a genome missing many more than the others: check its assembly or annotation |
| Singletons per genome | `_singletons` | pangenes found in only one genome | a genome with many more may be divergent, or carry contamination or annotation problems |
| POCS heatmap | `_POCS_heatmap` | percentage of conserved sequences between every pair of genomes, from `get_pangenes` (preliminary set only; step 07 does not change it) | how related the genomes are |
| InterPro enrichment | `_domain_enrichment` | InterPro entries found more often in core or cloud than in the other pangenes (one-sided Fisher test, Benjamini-Hochberg). Dot size: share of the category's pangenes with the entry; colour: adjusted p-value | an entry dominating the cloud can be an artefact (a transposable element, an organelle gene, a gene annotated in pieces) |

### Interactive explorer

An HTML page to browse every pangene: its occupancy and category, a strip with one cell per genome, and, for each gene,
its IDs, panel gene, model origin, ORF status, Liftoff coverage and identity, and protein domains. A locus view draws
every gene within 10 kb of the pangene, in each strain that has it (the same window in all of them, so the gene context
can be compared). The strains without it are seen in the pangenome browser.

- **Gene quality** comes from step 05's `orf_status`: High (valid, or repaired by miniprot), Low (invalid ORF or
  pseudogene), Not judged (ORF broken by the assembly: N in the CDS or contig end).
- **PSAURON**, once per pangene: the score of its panel gene's model (from `winners.tsv`). Below 0.5 is shown in red.
- **Invalid ORF in X of Y strains**, once per pangene: a strain counts as invalid when all its genes in the pangene are
  Low.
- **Absences**: every missing strain shows why the gene is missing (`absences_explained.tsv` of the safe fix, or of the
  removal for the final set).
- **Search**: pangene ID, gene ID, description, or protein domain (InterPro ID or name, Pfam ID): e.g. `IPR007226` or
  `SRS domain` lists every pangene with a gene carrying that domain.
- **Genes in no pangene** (left out by `get_pangenes`, or excluded from the pangenome in part 4), genes on contigs with
  no pangene included: grouped by panel gene, each group is an entry of its own, marked "in no pangene", found by the
  search and opened like a pangene (all strains, its genes framed). They count nowhere.
- **Filters**: category, occupancy, quality, strains present or absent, two genes in one strain, PSAURON < 0.5, and
  (preliminary page only) low confidence (invalid ORF in more than 50% of the strains, reference pseudogenes excluded).
- **Pangenome browser** (button at the top, or "Browser" in a pangene's panel): walk along the contigs of a guide strain
  and see every strain at once. Each column is one gene of the guide strain, in its order, with the pangene's category
  on top; each row is a strain: a filled cell when it has the gene (dark = valid ORF, light = invalid, "N" = a gap in
  the CDS, yellow = cut at a contig end), an empty cell
  when it does not (hover for the reason; "?" = maybe an assembly problem), and a red mark where the strain's gene is
  not next to the previous one (another contig, or genes in between: a rearrangement or an assembly break). Between two
  guide genes, **inserted columns** (italic orange ID) show the pangenes the guide strain does not have but other strains
  carry there, so every pangene of the region is seen; the guide strain only sets the order. Move with
  the arrows, the slider or the keyboard; type a pangene or gene ID to jump to it; click a column to open its pangene.
  **Highlight domain**: type an InterPro/Pfam ID or name; the columns whose pangene has it get an orange frame, the
  count of such pangenes is shown, and "match" buttons jump to the previous/next gene of the guide strain with it. With
  **show only these**, the browser shows only the pangenes with that domain, side by side: first those of the guide
  strain, in its order, then the ones it does not have.
- Genes and pangenes can be **flagged** on the page; the flags stay in the browser and can be exported.

Three optional layers are added when their input is there:

- **Original models** (step 05's gene fate and `input/native/`): under each strain, its own annotation, so a model that
  step 05 changed, joined or added can be compared with the original. An original gene with no final gene is drawn in
  red ("removed"), with the reason when hovered: why step 05 left it out, or, in the final page, the list of part 4
  that removed it (e.g. curation, TE). Built by [`build_source_layer.py`](scripts/build_source_layer.py).
- **Transposable elements** (`TE_DIR`: step 02's RepeatMasker output, one `RM_<genome>/` folder each): a line under each
  strain with the TEs of the window, named. Built by [`build_te_layer.py`](scripts/build_te_layer.py).
- **Orthogroups** (`OG_TABLE`: VEuPathDB's "Genes by Ortholog Pattern" table of the reference): the OrthoMCL group of
  each pangene's reference gene, and the pangenes that share one. Added by [`add_orthogroups.py`](scripts/add_orthogroups.py).
  Only the reference's genes are looked up (by their gene IDs), so a pangene with no reference gene gets no
  orthogroup. **To do:** give every gene of every strain its OrthoMCL orthogroup (e.g. by assigning each protein to
  the VEuPathDB OrthoMCL groups), so that every pangene can have one; the scripts for this step are not written yet.

**Strain names (optional).** When the strains are named by their assembly accession (e.g. `GCA_040126065.1`), a table
`07_Post_GET_PANGENES_adjusting/strain_names.tsv` (old name, new name, tab-separated) shows them under their strain
names in the explorer and the figures, and part 4 writes a copy of the final matrix with the new names
(`pangene_matrix_genes_strain_names.tr.tab`). [`rename_strains.py`](scripts/rename_strains.py) does the renaming. The
pipeline's own folders, GFFs and tables keep the original names. Without the table nothing changes.

```
old	new
GCA_040126065.1	AF100-9B
```

The **preliminary** page is for curation: it has the gene quality (High / Low), the flags and the PSAURON filter. The
**final** page is for users: no quality classes and no flags; each gene shows whether its ORF is valid; the description
and the protein domains are shown once, at the top of the pangene (the domains come from the protein of the source
models, not from each gene's own protein); and the genes left out of the pangenome (part 4) are drawn in red in the
locus view, with the reason.

The preliminary page (*T. gondii*): the list, with the pangene quality filter and the PSAURON filter; the top of a
pangene, with its quality (High) and the pangene flag; and the gene of one strain, with its quality and the curation
flag.

![Preliminary page: the list of pangenes](explorer_guide/images/preliminary_list.png)

![Preliminary page: the top of a pangene](explorer_guide/images/preliminary_pangene_top.png)

![Preliminary page: the gene of one strain](explorer_guide/images/preliminary_gene_card.png)

The final page is shown in [Part 5](#part-5--final-look).

[`build_explorer_job.sh`](scripts/build_explorer_job.sh) runs [`build_explorer_data.py`](scripts/build_explorer_data.py)
and [`build_explorer_locus_index.py`](scripts/build_explorer_locus_index.py), then packs
[`explorer.html`](scripts/explorer.html) and the data into one file, `explorer_standalone.html`. Copy that file to your
computer and open it with a double click (it is large, about 60 MB for *T. gondii*).

## Reference gene check

![A new gene on the reference check page](reference_check_guide/images/02_gene_supported.png)

A new ME49 gene on the page, with its intron junctions (purple arcs) and RNA-seq (green bars), and the choice made.
The [Reference gene check – user guide](reference_check_guide/README.md) explains the page. The page is a first help: the best place to
check a gene is still the genome browser (JBrowse) of VEuPathDB.

The reference is the genome most people use, so the genes the pipeline added to it (`not_in_reference_annotation=True`
in its step 05 GFF) are checked one by one. [`reference_check_job.sh`](scripts/reference_check_job.sh) writes one page,
`reference_check/reference_check.html` (open it in a browser). For each new gene it shows the gene with the reference
genes around it, its pangene and in how many strains, the strain its model came from, ORF, PSAURON, and the evidence:

- **RNA-seq** (`BW_ROOT`: a folder of bigWig files, searched with its subfolders; e.g. the VEuPathDB RNA-seq
  datasets): the coverage over the gene's CDS on the gene's own strand (stranded samples only: first strand = +,
  second strand = -; the reads of a gene on the other strand do not count), summed over the samples, as a percentile of 500 of the reference's own
  genes, and its profile (before the gene, 10 parts, after). The profile is drawn on one log scale for every gene, with
  a dashed line at the level of a typical reference gene (the median of the 500). Computed by
  [`region_rnaseq.py`](scripts/region_rnaseq.py); `RERUN_RNASEQ=no` reuses the table of an earlier run (to rebuild only
  the page);
- **intron junctions** (`JUNCTIONS`: e.g. VEuPathDB's `unifiedIntronJunction.gff`): drawn as arcs, with how many of the
  gene's introns a junction confirms;
- **BLAST link** (`GENOME_DIR`, `BLAST_EXCLUDE`): each new gene's protein, translated from `GENOME_DIR/<genome>.fna`,
  with a link that opens NCBI blastp against nr already filled in, the species of `BLAST_EXCLUDE` left out
  (e.g. `Aspergillus fumigatus (taxid:746128)`); the search runs only when BLAST is pressed on the NCBI page, so only
  the genes of interest are searched;
- **protein outside the species** (`PROTEIN_EVIDENCE`): the best blastp hit of the gene's protein in NCBI nr, the
  reference's own species left out, or "no hit" (e-value 1e-5); the page can show only the genes with or without a hit.
  The table comes from [`protein_evidence_job.sh`](scripts/protein_evidence_job.sh), run before the page:
  [`protein_evidence.py`](scripts/protein_evidence.py) searches NCBI over the internet (`blastp -remote`, no local
  database), 10 proteins at a time, and writes only `reference_check/protein_evidence/`. It can take many hours; a
  job that stops is submitted again and goes on from the batches already done.

All are optional: without `BW_ROOT`, `JUNCTIONS`, `GENOME_DIR` or `PROTEIN_EVIDENCE` the page is built without that
evidence. `protein_evidence_job.sh` needs internet access from the compute node.
`OUT_PAGE` gives the page another name, so an earlier page is kept.

For each gene one choice is made on the page: keep it as a new gene, delete it everywhere, leave it out of the
pangenome, or delete it in the reference only, with a reason picked from a list (supported, no support, mitochondrial
gene, model to fix, extra copy / gene family, overlaps a reference gene, other) and an optional detail. A gene with no
choice is deleted like a false gene (its whole panel gene, from every strain). The page says what each choice does to
the gene. The choices stay in the browser until the table is downloaded (`curation_actions.tsv`, every gene in it). Save it as `07_Post_GET_PANGENES_adjusting/curation_actions.tsv`: part 4 reads it together with
the lists. Or write the decisions in the lists by hand and run part 4 with `ACTIONS=`, so that no table is read
(e.g. when the decisions changed after the page). `delete` (A): the whole panel gene is removed from every strain, original genes included.
`exclude_from_pangenome` (B): every gene of the panel gene leaves the pangenome and its counts, but stays in the GFFs,
tagged `pangenome_excluded=<reason>` (`pangenome_final/excluded_from_pangenome.tsv`). `delete_in_genome` (C): the gene is not kept in its reference
annotation but stays in the pangenome and in these GFFs, tagged `reference_excluded=<reason>`
(`pangenome_final/reference_excluded.tsv`, for the official GFF). `keep_new_gene`: nothing.

## Part 4 – Removal after curation

Some genes are found to be false only by looking at the pangenome (the lists of part 2, the figures and the explorer of
part 3, the reference gene check): a transposable element that earlier clean-ups missed, a gene with an invalid ORF in
almost every strain and a low PSAURON score. Sometimes the gene is false in every strain; sometimes only one strain's
gene is wrong (e.g. a merged model where the strain has two good genes). Other genes are real but must not be counted
in the pangenome (e.g. pieces of the mitochondrial genome inserted in the nuclear genome). So there are three lists,
all in `07_Post_GET_PANGENES_adjusting/`:

- `genes_to_remove.tsv`: whole panel genes, removed from every strain (every gene with that `basket=`):

  ```
  panel_gene	reason
  Afu1g00140	Pogo/Tc1-mariner transposable element (HTH CenpB + HTH Psq domains)
  ```

- `strain_genes_to_remove.tsv`: single genes, removed only in their strain; the other strains keep theirs:

  ```
  gene	reason
  TGGT1_221865	example: false gene in GT1 only
  ```

  A gene of the final GFFs has two IDs, e.g. `ID=Tg7g00097;original_ID=TGGT1_221865`: the new ID given by step 05,
  and the strain's own gene ID. Either can be used. The strain's own ID is safer: the new IDs are numbered along each
  genome, so if step 05 is ever run again (e.g. with new rules or a new annotation) they can change, and a list of new
  IDs would point to other genes. Genes that step 05 added to a strain (taken from another strain) have no
  `original_ID`; for them, use the new ID, and check the list again after any new run of step 05.

- `genes_to_exclude_from_pangenome.tsv`: whole panel genes left out of the pangenome only. Their genes leave the matrix
  and every count, but stay in the GFFs, tagged `pangenome_excluded=<reason>`, so the explorer still shows them:

  ```
  panel_gene	reason
  TGRH88_085440	mito: cytochrome b
  ```

A list with only its header changes nothing.

[`remove_panel_genes.py`](scripts/remove_panel_genes.py) removes the genes with all their features (mRNA, exons, CDS)
from a copy of the GFFs and from the matrix, and takes the genes of the exclusion list out of the matrix only; a
pangene left empty disappears. An ID of the single-gene list found in no GFF is listed in `genes_not_found.tsv`. Empty
lists (header only) copy everything unchanged, so `pangenome_final/` always holds the result of step 07.
`removed_genes.tsv` says which list removed each gene; `excluded_from_pangenome.tsv` lists every gene left out of the
pangenome, with its model origin and reason.

**Pangene IDs.** `get_pangenes` names each pangene after one of its genes (e.g. `Tg8g00699`), which reads as if the
pangene belonged to that strain, and changes when step 06 is run again. The final pangenes are renamed
`<PANGENE_PREFIX>_pan_00001`, `_00002`... in the order of the matrix (`get_pangenes`' reference genome, along its
contigs). `pangene_ids.tsv` gives the old name of each, its number of strains and genes, and its panel genes. The final
matrix, `absences_explained.tsv` and the tables below use the new IDs.

Then [`final_tables.py`](scripts/final_tables.py) describes the final pangenome:

- `gene_status.tsv`: every gene of every strain in step 05's GFFs, with its pangene in `get_pangenes` and in the final
  pangenome, and its status: `in_pangenome`, `moved` (by part 2), `not_in_pangenome` (CDS < 100 bp),
  `removed_curation` (whole panel gene), `removed_curation_gene` (single gene) or `excluded_from_pangenome`, with the
  reason. The panel genes that step 05 itself removed or left out are listed too, so every
  gene of the panel is accounted for.
- `singletons.tsv`: every pangene present in one strain, and why it is alone: `tandem_copy`, `split_by_get_pangenes`
  (its panel gene also has genes in other pangenes; the safe fix could not join them), `inside_another_gene`,
  `partial` (CDS < 50% of its panel gene), `unique_to_strain` (no other strain has its panel gene) or `other`.

## Part 5 – Final look

The same figures and explorer as in part 3, on `pangenome_final/` (`SET=final`, the default); the final explorer is the
version for users (see [Interactive explorer](#interactive-explorer)). Repeat parts 4 and 5 after each round of curation.

The final page (*T. gondii*): the same pangene, now `Tg_pan_01197`, with its description and domains at the top, no
quality classes and no flags, and its locus view. The [Pangenome Explorer – user guide](explorer_guide/README.md) describes every part of
the final page, with more screenshots.

![Final page: the top of a pangene](explorer_guide/images/04_pangene_top.png)

![Final page: the locus view](explorer_guide/images/05_locus.png)

`pangenome_final/gff/` holds working GFFs: the step 05 models, with the tags of steps 05 and 07. The final annotation is made from
them in [step 08](../08_Final_GFF/README.md).

## Running it

From the project folder on the cluster, in this order. Each log folder must exist before submitting.

```bash
# 1. inspection
mkdir -p pangenome_inspection/logs
sbatch 07_Post_GET_PANGENES_adjusting/scripts/inspect_pangenome_job.sh
# 2. safe fix (after 1)
mkdir -p pangenome_fixed/logs
sbatch 07_Post_GET_PANGENES_adjusting/scripts/fix_pangenome_safe_job.sh
# 3. preliminary look (after 2)
mkdir -p pangenome_plots/logs pangenome_explorer/logs
sbatch --export=ALL,SET=preliminary 07_Post_GET_PANGENES_adjusting/scripts/plot_all_job.sh
sbatch --export=ALL,SET=preliminary 07_Post_GET_PANGENES_adjusting/scripts/build_explorer_job.sh
# reference gene check (after 2): the RNA-seq bigWig folder and the intron junctions GFF are optional
mkdir -p reference_check/logs
sbatch --export=ALL,BW_ROOT=/path/to/bigwigs,JUNCTIONS=unifiedIntronJunction.gff 07_Post_GET_PANGENES_adjusting/scripts/reference_check_job.sh
#    the step 05 and part 2 folders of the same run, when they are not basket_annotation/ and pangenome_fixed/:
#    sbatch --export=ALL,STEP05_DIR=<step 05 folder>,FIXED_DIR=<part 2 folder> ...
#    only the page again (the RNA-seq table of the previous run is reused):
#    sbatch --export=ALL,RERUN_RNASEQ=no ...
#    protein evidence (optional, before the page; EXCLUDE_TAXID = the reference's species, e.g. 746128 A. fumigatus,
#    5811 T. gondii; submit again to go on if it stops):
#    sbatch --export=ALL,EXCLUDE_TAXID=<taxid> 07_Post_GET_PANGENES_adjusting/scripts/protein_evidence_job.sh
#    then the page with it, under another name:
#    sbatch --export=ALL,RERUN_RNASEQ=no,PROTEIN_EVIDENCE=reference_check/protein_evidence/protein_evidence.tsv,OUT_PAGE=reference_check_protein.html ...
# -> curation: add the genes to remove to 07_Post_GET_PANGENES_adjusting/genes_to_remove.tsv (whole panel genes)
#    or strain_genes_to_remove.tsv (single genes); the genes to leave out of the pangenome only to
#    genes_to_exclude_from_pangenome.tsv
# 4. removal after curation (after 2; again after each change to the lists). PANGENE_PREFIX names the final pangenes
#    (e.g. Tg -> Tg_pan_00001); ACTIONS= leaves out a curation_actions.tsv saved in the folder
mkdir -p pangenome_final/logs
sbatch --export=ALL,PANGENE_PREFIX=Tg,ACTIONS= 07_Post_GET_PANGENES_adjusting/scripts/remove_panel_genes_job.sh
# 5. final look (after 4)
sbatch 07_Post_GET_PANGENES_adjusting/scripts/plot_all_job.sh
sbatch 07_Post_GET_PANGENES_adjusting/scripts/build_explorer_job.sh
```

The settings are at the top of each job, between `===== Edit these =====` lines, and can also be given on the command
line, e.g. `sbatch --export=ALL,LABEL=Afum ...`. The defaults expect step 05's outputs in `basket_annotation/` and
`get_pangenes`' output in `for_get_pangenes_basket_prefer_panel_model_pangenes/`.

## Outputs

```
pangenome_plots/{preliminary,final}/      <LABEL>_<SET>_*.png (+ .tsv/.log with the numbers)
pangenome_explorer/{preliminary,final}/   explorer_standalone.html (+ explorer.html, data/)
pangenome_inspection/
├── absence_audit/      absences.tsv, summary.tsv
└── multigene/          multigene_pangenes.tsv, summary.tsv
pangenome_fixed/        pangene_matrix_genes.tr.tab/.tab, moves.tsv, skipped_moves.tsv, joins.tsv,
                        basket_splits_not_joined.tsv, absences_explained.tsv, summary.tsv,
                        pangene_lists/ (pangene_quality.tsv, genes_not_in_pangenome.tsv, summary.tsv)
reference_check/        reference_check.html, <reference>_region_rnaseq.tsv, protein_evidence/ (protein_evidence.tsv)
pangenome_final/        gff/ (working GFFs), pangene_matrix_genes.tr.tab/.tab, pangene_ids.tsv, removed_genes.tsv,
                        pangene_matrix_genes_strain_names.tr.tab (with strain_names.tsv),
                        excluded_from_pangenome.tsv, genes_not_found.tsv, absences_explained.tsv, gene_status.tsv,
                        singletons.tsv, summary.tsv, summary_tables.tsv
```

`pangenome_final/pangene_matrix_genes.tr.tab` is the final pangenome. `pangenome_final/gff/` is the input of
[step 08](../08_Final_GFF/README.md), which makes the final annotation.

## Scripts

| Part | Job | Scripts |
|---|---|---|
| 3, 5 | [`plot_all_job.sh`](scripts/plot_all_job.sh) | [`plot_all.sh`](scripts/plot_all.sh), [`plot_pancore.py`](scripts/plot_pancore.py), [`plot_occupancy.py`](scripts/plot_occupancy.py), [`plot_softcore_missing.py`](scripts/plot_softcore_missing.py), [`plot_singletons.py`](scripts/plot_singletons.py), [`plot_pocs.py`](scripts/plot_pocs.py), [`plot_domain_enrichment.py`](scripts/plot_domain_enrichment.py), [`rename_strains.py`](scripts/rename_strains.py) |
| 3, 5 | [`build_explorer_job.sh`](scripts/build_explorer_job.sh) | [`build_explorer_data.py`](scripts/build_explorer_data.py), [`build_explorer_locus_index.py`](scripts/build_explorer_locus_index.py), [`explorer.html`](scripts/explorer.html), [`build_source_layer.py`](scripts/build_source_layer.py), [`build_te_layer.py`](scripts/build_te_layer.py), [`add_orthogroups.py`](scripts/add_orthogroups.py), [`rename_strains.py`](scripts/rename_strains.py), [`make_standalone_explorer.py`](scripts/make_standalone_explorer.py), [`serve_explorer.py`](scripts/serve_explorer.py) (to open `explorer.html` with its data folder instead of the single file) |
| 1 | [`inspect_pangenome_job.sh`](scripts/inspect_pangenome_job.sh) | [`audit_absences.py`](scripts/audit_absences.py), [`inspect_multigene_pangenes.py`](scripts/inspect_multigene_pangenes.py) |
| 2 | [`fix_pangenome_safe_job.sh`](scripts/fix_pangenome_safe_job.sh) | [`fix_pangenome_safe.py`](scripts/fix_pangenome_safe.py), [`summarize_pangenes.py`](scripts/summarize_pangenes.py) |
| reference check | [`reference_check_job.sh`](scripts/reference_check_job.sh) | [`region_rnaseq.py`](scripts/region_rnaseq.py), [`build_reference_check.py`](scripts/build_reference_check.py), [`reference_check_template.html`](scripts/reference_check_template.html) |
| protein evidence (optional) | [`protein_evidence_job.sh`](scripts/protein_evidence_job.sh) | [`protein_evidence.py`](scripts/protein_evidence.py) |
| 4 | [`remove_panel_genes_job.sh`](scripts/remove_panel_genes_job.sh) | [`remove_panel_genes.py`](scripts/remove_panel_genes.py), [`final_tables.py`](scripts/final_tables.py), [`rename_strains.py`](scripts/rename_strains.py) |

[↑ Back to top](#07--post-get_pangenes-adjusting-and-visualisation) &nbsp;|&nbsp; [Back to Pangenome Pipeline](../README.md)
