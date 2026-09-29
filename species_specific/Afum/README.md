# *Aspergillus fumigatus*-specific scripts

One-off scripts written for the *A. fumigatus* run. They are not general pipeline steps. The folders mirror the pipeline's step layout, and each script is submitted from the project folder like the others. See the [*A. fumigatus* report](../../results/Afum/README.md) for when and why each one was used.

| Script | What it does |
|---|---|
| [`step02_te_annotation_clean_up/fix_repeatmasker_after_contig_drop.slurm.sh`](step02_te_annotation_clean_up/fix_repeatmasker_after_contig_drop.slurm.sh) · [`.py`](step02_te_annotation_clean_up/fix_repeatmasker_after_contig_drop.py) | Brings existing RepeatMasker output back in line after contigs were dropped in step 01, instead of re-running RepeatMasker. |
