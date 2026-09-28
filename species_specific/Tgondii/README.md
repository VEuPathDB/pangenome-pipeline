# *Toxoplasma gondii*-specific scripts

One-off scripts written for the *T. gondii* run. They are not general pipeline steps. The folders mirror the pipeline's step layout, and each script is submitted from the project folder like the others. See the [*T. gondii* report](../../results/Tgondii/README.md) for when and why each one was used.

| Script | What it does |
|---|---|
| [`01_Genome_Selection_QC/scripts/drop_short_contigs.slurm.sh`](01_Genome_Selection_QC/scripts/drop_short_contigs.slurm.sh) | Drops every contig below `MIN_LEN` from `genomes/` and `gff_originals/`, using the general [`drop_contigs.py`](../../src/pangenome_pipeline/step01_genome_selection_qc/drop_contigs.py). It backs up every affected file first. |
| [`09_Redundant_Contig_Removal/scripts/make_removal_list.py`](09_Redundant_Contig_Removal/scripts/make_removal_list.py) | Builds the list of small, gene-carrying redundant contigs to remove after the redundancy check. |
| [`09_Redundant_Contig_Removal/scripts/remove_redundant_contig_genes.py`](09_Redundant_Contig_Removal/scripts/remove_redundant_contig_genes.py) | Removes those contigs' genes from the GFFs and the get_pangenes matrices, writing to new folders. |

**Not in this repo yet:** these scripts are referenced but were never committed. They are still in the *T. gondii* species project:

- `01_Genome_Selection_QC/scripts/summarize_dropped_genes.py`, called by `drop_short_contigs.slurm.sh`
- `contig_redundancy_final.slurm.sh`
- `09_Redundant_Contig_Removal/scripts/remove_redundant_contig_genes.slurm.sh`
- `09_Redundant_Contig_Removal/scripts/plot_nonredundant_job.sh`
