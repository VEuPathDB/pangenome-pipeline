"""Python code of the pangenome pipeline.

One subpackage per numbered step (step01_genome_selection_qc, ..., step08_pangenome_visualisation)
plus common/ (helpers shared across steps) and alternative_get_panloci/. The step folders at the repo
root keep each step's README and SLURM job scripts; the jobs run the code here with
`python3 -m pangenome_pipeline.<step>.<module>`, with src/ put on PYTHONPATH by lib/slurm_common.sh.
"""
