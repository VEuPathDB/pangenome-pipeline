# *Toxoplasma gondii*-specific scripts

One-off scripts written for the *T. gondii* run. They are not general pipeline steps. The folders mirror the pipeline's step layout, and each script is submitted from the project folder like the others. See the [*T. gondii* report](../../results/Tgondii/README.md) for when and why each one was used.

| Script | What it does |
|---|---|
| [`step01_genome_selection_qc/drop_short_contigs.slurm.sh`](step01_genome_selection_qc/drop_short_contigs.slurm.sh) | Drops every contig below `MIN_LEN` from `genomes/` and `gff_originals/`, using the general [`drop_contigs.py`](../../src/pangenome_pipeline/step01_genome_selection_qc/drop_contigs.py). It backs up every affected file first. |
| [`step09_redundant_contig_removal/make_removal_list.py`](step09_redundant_contig_removal/make_removal_list.py) | Builds the list of small, gene-carrying redundant contigs to remove after the redundancy check. |
| [`step09_redundant_contig_removal/remove_redundant_contig_genes.py`](step09_redundant_contig_removal/remove_redundant_contig_genes.py) | Removes those contigs' genes from the GFFs and the get_pangenes matrices, writing to new folders. |

**Not in this repo yet:** these scripts are referenced but were never committed. They are still in the *T. gondii* species project:

- `step01_genome_selection_qc/summarize_dropped_genes.py`, called by `drop_short_contigs.slurm.sh` (it looks for it there; set `SUMMARIZE_PY` to use a copy elsewhere - without it the job skips the gene-cost summary)
- `contig_redundancy_final.slurm.sh`
- `step09_redundant_contig_removal/remove_redundant_contig_genes.slurm.sh`
- `step09_redundant_contig_removal/plot_nonredundant_job.sh`
