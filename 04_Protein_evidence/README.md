# 04 – Protein evidence (InterProScan and PSAURON)

***Purpose:*** Make sure every genome carries InterPro/PFAM protein-domain annotation in its
own original GFF, so that domain evidence can be used consistently to validate gene models
later in the pipeline, and score every protein with PSAURON (used in step 05).

## Running InterProScan across the panel

A 4-step chain: protein extraction, then a chunk manifest (splits each genome's proteins into
`CHUNK_SIZE`-sized chunks so InterProScan runs on manageable pieces instead of one huge job per
genome), then the array job (runs InterProScan itself, one task per chunk), then the merge.

**Which GFF to extract proteins from:** default is **`gff_cleaned/`** — the GFFs step 02's TE
kill list already filtered. For a species with **no TE library** (step 02 skipped
entirely), use **`gff_originals/`** instead — by step 04 it has already been edited in place
by step 01's organelle + redundant-contig removal, so it's the right fallback, not the raw
downloaded GFF. Set this via `GFF_INDIR` at the top of `extract_proteins.slurm.sh` (script 1
below).

**Dependencies:**
- **AGAT** (steps 1 and 5 below) — hardcoded here with an env called `bioinfo` (adapt to your
  own).
- **InterProScan + Java** (step 3) — via `module load`
  (`interproscan/5.66-98.0`, `adoptopenjdk/11.0.12+7`), not micromamba.
- Steps 2 and 4 (manifest + merge) only need plain `python3`, already on the PATH - no env to
  set up for those.

Scripts, in the order they run:

1. [`scripts/extract_proteins.slurm.sh`](scripts/extract_proteins.slurm.sh) — **new.** Runs
   AGAT (`agat_sp_extract_sequences.pl -p -t cds`) on every GFF in `GFF_INDIR` against its
   matching FASTA in `GENOMES`, writing one translated, spliced protein FASTA per genome to
   `OUTDIR`.
   - `GFF_INDIR="gff_cleaned"` — step 02's TE-kill-list output. Point at `gff_originals` only
     if you're deliberately skipping step 02.
   - `GENOMES="genomes"` — the step-01 assemblies folder.
   - `OUTDIR="interpro/proteins"` — `<genome_label>.faa` land here.
   - `### ADAPT ###` block — put `agat_sp_extract_sequences.pl` on the PATH (the `bioinfo` env
     used for AGAT elsewhere in this project).
   ```bash
   sbatch 04_Protein_evidence/scripts/extract_proteins.slurm.sh
   ```

2. [`scripts/build_interpro_manifest.slurm.sh`](scripts/build_interpro_manifest.slurm.sh) —
   splits each genome's protein FASTA (`.faa`, from step 1 above — the genome assembly itself
   stays a `.fna` in `GENOMES` the whole time) into `CHUNK_SIZE`-protein chunks and writes one
   manifest row per (genome, chunk): `task_id, genome, chunk_num, start_seq, end_seq`. The
   array job below reads this to know what each task should process.
   - `PROTEIN_DIR="interpro/proteins"`, `ORDER="genome_order.txt"` (this project's canonical
     genome list, same one used in steps 01/02/04), `CHUNK_SIZE=1000`,
     `MANIFEST_OUT="interpro/interpro_manifest.tsv"`.
   ```bash
   sbatch 04_Protein_evidence/scripts/build_interpro_manifest.slurm.sh
   ```
   Prints `Total chunks: N` at the end — that's the array size for the next step.

3. [`scripts/interproscan_array_job.slurm.sh`](scripts/interproscan_array_job.slurm.sh) — one
   SLURM array task per manifest row; each task extracts its own chunk's sequences from the
   genome's `.faa` and runs `interproscan.sh -f tsv -goterms -iprlookup -pa -dp` on them.
   - `-dp` (disable precalculated match lookup) is required because compute nodes on this
     cluster can't reach EBI — every sequence is computed locally instead of checked against
     EBI's online cache, which is slower and already accounted for in the `--time` budget.
   - `### ADAPT ###` block — `module load` (or equivalent) for Java + InterProScan.
   - Output: `interpro/interproscan_results/<genome>/chunk_<N>_interproscan.tsv`.
   - **Never hand-edit the `#SBATCH --array=1-1%20` line** in this script to the real chunk
     count — that line is only a smoke-test-safe default. `--array=...` passed on the `sbatch`
     command line always overrides it, so compute the real count from the manifest at submit
     time instead — it can never go stale that way:
   ```bash
   mkdir -p interpro/logs   # SLURM needs this to exist at job start
   N=$(($(wc -l < interpro/interpro_manifest.tsv) - 1))   # -1 for the manifest's header row
   sbatch --array=1-${N}%20 04_Protein_evidence/scripts/interproscan_array_job.slurm.sh
   ```
   (`%20` caps concurrent tasks — adjust to your cluster's fair-use policy.) To smoke-test a
   single task first instead:
   ```bash
   sbatch --array=1-1 04_Protein_evidence/scripts/interproscan_array_job.slurm.sh
   ```

4. [`scripts/merge_interpro_job.slurm.sh`](scripts/merge_interpro_job.slurm.sh) (calls
   [`scripts/merge_interpro_results.py`](scripts/merge_interpro_results.py)) — run only once
   **every** array task has finished. For each genome, collapses all its
   `chunk_*_interproscan.tsv` files into one `<genome>_interpro.tsv`:
   `protein_id, interpro_ids, pfam_ids, n_domain_hits` (semicolon-separated IDs when a protein
   has more than one hit) — this is the lookup table step 5 below reads from.
   - `RESULTS_DIR="interpro/interproscan_results"`, `PROTEIN_DIR="interpro/proteins"` (used
     only for the completeness check below), `ORDER="genome_order.txt"`,
     `OUT_DIR="interpro/interpro_lookup"`.
   - Two completeness checks per genome, printed to the log: a chunk number with no matching
     TSV file (`WARNING` — InterProScan never produced output for it, a real problem, check
     that chunk's own job log) and a protein present in the genome's `.faa` but absent from
     every chunk TSV (`NOTE`, not `WARNING` — InterProScan's `-f tsv` format only ever writes a
     row when an analysis finds a hit, so this is almost always just a protein with zero domain
     hits across every analysis, not a failure; there's no way to tell that apart from a
     genuinely malformed sequence without spot-checking one by hand). Only the first (missing
     chunk file) means the lookup table is actually incomplete.
   ```bash
   sbatch 04_Protein_evidence/scripts/merge_interpro_job.slurm.sh
   ```

5. [`scripts/update_interpro_in_gff.slurm.sh`](scripts/update_interpro_in_gff.slurm.sh) (calls
   [`scripts/update_interpro_in_gff.py`](scripts/update_interpro_in_gff.py)) — writes the
   `Dbxref=InterPro:.../PFAM:...` evidence into each genome's own
   GFF — the field step 05's domain-evidence check actually reads. One genome at a time: for
   each CDS line, strips any `InterPro:`/`PFAM:` entries the original annotation already had
   (could be from an older InterProScan run) and replaces them with this run's own result,
   matched by the CDS's `Parent=` — the mRNA ID, which is also step 4's lookup key, since
   that's the same ID AGAT used as the FASTA header when it extracted that protein in step 1.
   Everything else in `Dbxref=` (e.g. `NCBI_GP:...`) is left alone.
   - `GFF_INDIR="gff_cleaned"`, `LOOKUP_DIR="interpro/interpro_lookup"`,
     `ORDER="genome_order.txt"`, `GFF_OUTDIR="gff_with_interpro"` — new GFFs land here,
     `GFF_INDIR` itself is never modified.
   ```bash
   sbatch 04_Protein_evidence/scripts/update_interpro_in_gff.slurm.sh
   ```
   Log line per genome: how many CDS lines had old evidence stripped, how many got fresh
   evidence, and the match percentage — check this before trusting the output for a genome.

   **Before or after Liftoff.** This pipeline runs InterProScan **before** Liftoff, on each genome's
   own original annotation (this step). The alternative is to run it **after** Liftoff, on the
   transferred/final annotation: the original scripts work that way.

## Scoring every protein with PSAURON

PSAURON scores how much a sequence looks like a real protein-coding gene (0 to 1). Step 05 uses the scores to decide
which genes not annotated in the reference can enter the panel, and to settle overlaps (`PSAURON_DIR="psauron"` in its
job script).

**Dependencies:** a Singularity container with PSAURON (`.sif`; we used PSAURON 1.0.8), run with `singularity exec`.
Set its path in the script's `### ADAPT ###` block.

[`scripts/psauron_array_job.slurm.sh`](scripts/psauron_array_job.slurm.sh) runs PSAURON on every protein extracted by
[`scripts/extract_proteins.slurm.sh`](scripts/extract_proteins.slurm.sh) (script 1 above), one array task per genome.
It only needs those proteins, so it can run at the same time as InterProScan.
- `PROTEIN_DIR="interpro/proteins"`, `ORDER="genome_order.txt"`, `OUT_DIR="psauron"`; `### ADAPT ###` — the container
  path (`singularity` on the PATH).
- Output: `psauron/<genome>.csv` (`description,psauron_is_protein,in-frame_score`; `description` is the mRNA ID).
  The log gives the number of scores per genome: it should equal the number of proteins.
- A genome whose CSV is already complete is skipped, so the same command resumes a run.
```bash
mkdir -p psauron/logs
N=$(grep -v '^#' genome_order.txt | awk 'NF' | wc -l)
sbatch --array=1-${N} 04_Protein_evidence/scripts/psauron_array_job.slurm.sh
```

[Back to top level](../README.md)
