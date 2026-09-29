# Moving an existing project to the new layout

Until September 2026 the pipeline was cloned **into** each species project, every step had a
numbered folder (`01_Genome_Selection_QC/scripts/`, ...), and jobs wrote their results into the
project root or into those step folders. Now the pipeline checkout and the project are separate
folders, the job scripts are in [`workflow/`](../workflow/), the step pages in `docs/`, and every
job writes under `$OUTPUT_DIR/stepNN_<name>/` (by default `output/` in the project folder) - see
[getting the pipeline into a project](../README.md#getting-the-pipeline-into-a-project).

A project that already has results (Afum, Tgondii) keeps working once its results are moved to where
the jobs now look for them. This page lists every path that changed and gives a script that does the
move.

## 1. Get the new pipeline

Either clone the pipeline somewhere new and point `PIPELINE_DIR` at it:

```bash
git clone https://github.com/adriludwig/pangenome-pipeline.git ~/pangenome-pipeline
echo 'export PIPELINE_DIR=~/pangenome-pipeline' >> ~/.bashrc
cp <project>/pipeline_config.sh ~/pangenome-pipeline/pipeline_config.sh   # your cluster settings
```

or `git pull` inside the clone that sits in the project folder. Git then moves the job scripts and
READMEs out of the numbered folders, but it leaves the **untracked results** in them (for example
`02_TE_annotation_Clean_Up/final_TE_candidates_after_MCHelper/`), so run the script below either way.

## 2. Move the results

Run this from the project folder. It only moves what exists, never overwrites anything already at the
new place (it says so and leaves both), and prints every move:

```bash
mv_to() {   # mv_to <old path> <new path>
    local src="$1" dst="$2"
    [[ -e "$src" || -L "$src" ]] || return 0
    if [[ -e "$dst" ]]; then echo "SKIP  $src  ($dst already exists)"; return 0; fi
    mkdir -p "$(dirname "$dst")" && mv "$src" "$dst" && echo "MOVED $src -> $dst"
}
O=output
S1=$O/step01_genome_selection_qc  S2=$O/step02_te_annotation_clean_up  S3=$O/step03_denovo_genome_annotation
S4=$O/step04_interproscan  S5=$O/step05_cross_genome_annotation_integration  S6=$O/step06_get_pangenes
S7=$O/step07_post_get_pangenes_adjusting  S8=$O/step08_pangenome_visualisation

# step 01
for d in species_check synteny contig_size_profile organelle_screen contig_redundancy; do
    mv_to "$d" "$S1/$d"; mv_to "01_Genome_Selection_QC/$d" "$S1/$d"
done
for f in drop_contigs.txt drop_short_contigs.txt dropped_genes_summary.tsv backup_before_short_contig_drop; do
    mv_to "01_Genome_Selection_QC/$f" "$S1/$f"
done
# step 02 - RepeatMasker folders, the TE-library pipeline, the TE/gene overlap chain
for d in RM_* R_*; do mv_to "$d" "$S2/repeatmasker/$d"; done
for d in te_content_summary te_annotation_integration kill_lists gff_cleaned repeatmasker_contig_filtered te_aid_*; do
    mv_to "$d" "$S2/$d"
done
mv_to curated_TE_library.fa "$S2/curated_TE_library.fa"
for d in config repeatmodeler_more_genomes TE_candidate_screen_original_headers MCHelper_newTEs \
         final_TE_candidates_after_MCHelper TE_gap_screen; do
    mv_to "02_TE_annotation_Clean_Up/$d" "$S2/te_library/$d"
done
for f in 02_TE_annotation_Clean_Up/repbase_nt* 02_TE_annotation_Clean_Up/te_peptides*; do mv_to "$f" "$S2/te_library/$(basename "$f")"; done
for f in 02_TE_annotation_Clean_Up/*.fa 02_TE_annotation_Clean_Up/*.fasta; do mv_to "$f" "$S2/$(basename "$f")"; done
# step 03 - its whole working folder (inputs and outputs)
if [[ -d 03_denovo_Genome_Annotation ]]; then
    for f in 03_denovo_Genome_Annotation/*; do
        case "$(basename "$f")" in scripts|README.md) ;; *) mv_to "$f" "$S3/$(basename "$f")" ;; esac
    done
fi
# step 04
for d in proteins interproscan_results interpro_lookup; do mv_to "interpro/$d" "$S4/$d"; done
mv_to interpro/interpro_manifest.tsv "$S4/interpro_manifest.tsv"
mv_to gff_with_interpro "$S4/gff_with_interpro"
# step 05
mv_to kill_list_manual.tsv "$S5/kill_list_manual.tsv"
for d in annotation_check gff_curated gff_for_liftoff liftoff_output expanded_annotation \
         filtered_annotation tagged_annotation filter_qc_plots final_annotation_qc_plots; do
    mv_to "$d" "$S5/$d"
done
# steps 06-08
for d in for_get_pangenes_7a for_get_pangenes_7b for_get_pangenes_7a_pangenes for_get_pangenes_7b_pangenes; do
    mv_to "$d" "$S6/$d"
done
mv_to multigene_pangenes "$S7/multigene_pangenes"
mv_to pangenome_plots "$S8/pangenome_plots"
mkdir -p logs
```

Afterwards:

- **Re-run the two step 06 link jobs** (`link_for_get_pangenes_7a.slurm.sh` / `_7b`). Their outputs are
  symlinks with absolute paths into `filtered_annotation/`, `tagged_annotation/` and `genomes/`, so the
  ones that were moved now point at the old places. The link jobs take seconds and overwrite them.
- **Old SLURM logs** (`*_<jobid>.out` / `.err` in the project root, `<folder>/logs/` inside the moved
  folders) are not moved; new jobs log to `logs/`.
- Anything the script did not know about stays where it was; the numbered folders can be deleted once
  they hold nothing you want to keep.
- Folders you had pointed somewhere else by editing a job's settings are not covered - move those by
  hand, using the table below.

## 3. What moved where

Paths relative to the project folder; `output/stepNN_*` is `$OUTPUT_DIR/stepNN_*` (default `output/`).

**Pipeline files** (now in the pipeline checkout, not in the project):

| Old | New |
|---|---|
| `NN_<Step>/scripts/<job>` | `workflow/stepNN_<name>/<job>` |
| `02_TE_annotation_Clean_Up/slurm/<file>` | `workflow/step02_te_annotation_clean_up/te_library/<file>` |
| `Alternative_Get_PanLoci_pipeline/scripts/`, `scripts_no_slurm/` | `workflow/alternative_get_panloci/`, `workflow/alternative_get_panloci/no_slurm/` |
| `NN_<Step>/README.md` | `docs/stepNN_<name>.md` (`02_.../slurm/README.md` -> `docs/step02_te_library.md`) |
| `species_specific/<Sp>/NN_<Step>/scripts/<file>` | `species_specific/<Sp>/stepNN_<name>/<file>` |
| `pipeline_config.sh` in the project | `pipeline_config.sh` in the checkout; a project may keep its own for overrides |

**Results:**

| Step | Old | New |
|---|---|---|
| 01 | `species_check/`, `synteny/`, `contig_size_profile/`, `organelle_screen/`, `contig_redundancy/` | `output/step01_genome_selection_qc/<same>/` |
| 01 | `01_Genome_Selection_QC/drop_contigs.txt` (and Tgondii's `drop_short_contigs.txt`, `dropped_genes_summary.tsv`, `backup_before_short_contig_drop/`) | `output/step01_genome_selection_qc/<same>` |
| 02 | `RM_<genome>/` | `output/step02_te_annotation_clean_up/repeatmasker/RM_<genome>/` |
| 02 | `te_content_summary/`, `te_annotation_integration/`, `kill_lists/`, `gff_cleaned/`, `te_aid_*/`, Afum's `repeatmasker_contig_filtered/` | `output/step02_te_annotation_clean_up/<same>/` |
| 02 | `curated_TE_library.fa` | `output/step02_te_annotation_clean_up/curated_TE_library.fa` |
| 02 | `02_TE_annotation_Clean_Up/{config,repeatmodeler_more_genomes,TE_candidate_screen_original_headers,MCHelper_newTEs,final_TE_candidates_after_MCHelper}/`, `repbase_nt*`, `te_peptides*` | `output/step02_te_annotation_clean_up/te_library/<same>` |
| 03 | everything in `03_denovo_Genome_Annotation/` (inputs and outputs) | `output/step03_denovo_genome_annotation/` |
| 04 | `interpro/{proteins,interpro_manifest.tsv,interproscan_results,interpro_lookup}` | `output/step04_interproscan/<same>` |
| 04 | `gff_with_interpro/` | `output/step04_interproscan/gff_with_interpro/` |
| 05 | `kill_list_manual.tsv`, `annotation_check/`, `gff_curated/`, `gff_for_liftoff/`, `liftoff_output/`, `expanded_annotation/`, `filtered_annotation/`, `tagged_annotation/`, `filter_qc_plots/`, `final_annotation_qc_plots/` | `output/step05_cross_genome_annotation_integration/<same>` |
| 06 | `for_get_pangenes_7a/`, `_7b/`, `_7a_pangenes/`, `_7b_pangenes/` | `output/step06_get_pangenes/<same>/` |
| 07 | `multigene_pangenes/` | `output/step07_post_get_pangenes_adjusting/multigene_pangenes/` |
| 08 | `pangenome_plots/` | `output/step08_pangenome_visualisation/pangenome_plots/` |
| Get-PanLoci | the `RUN_DIR` you set | unchanged if you keep setting `RUN_DIR`; default is now `output/alternative_get_panloci/` |
| logs | `<job>_<id>.out` in the project, or `<output folder>/logs/` | `logs/<job>_<id>.out` |

**Renamed settings:** `OUTPUT_DIR` in `build_final_annotation_basic.sh` and
`check_final_annotations_job.sh` is now `FINAL_DIR` (`OUTPUT_DIR` is the pipeline-wide output
folder); `PROJECT_DIR` in the Get-PanLoci step 0 scripts is now `RUN_DIR`.

[Back to top level](../README.md)
