# *Aspergillus fumigatus*-specific scripts

One-off scripts written for the *A. fumigatus* run. They are not general pipeline steps. The folders mirror the pipeline's step layout, and each script is submitted from the project folder like the others. See the [*A. fumigatus* report](../../results/Afum/README.md) for when and why each one was used.

| Script | What it does |
|---|---|
| [`02_TE_annotation_Clean_Up/scripts/fix_repeatmasker_after_contig_drop.slurm.sh`](02_TE_annotation_Clean_Up/scripts/fix_repeatmasker_after_contig_drop.slurm.sh) · [`.py`](02_TE_annotation_Clean_Up/scripts/fix_repeatmasker_after_contig_drop.py) | Brings existing RepeatMasker output back in line after contigs were dropped in step 01, instead of re-running RepeatMasker. |
