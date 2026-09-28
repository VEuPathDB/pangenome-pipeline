#!/usr/bin/env python3
"""
Merges every chunk_N_interproscan.tsv (InterProScan's own '-f tsv' output format) for one
genome into a single per-protein domain lookup table.

InterProScan TSV columns (no header row in the file itself):
  1  Protein accession
  2  Sequence MD5 digest
  3  Sequence length
  4  Analysis (source database, e.g. Pfam, PANTHER, ProSiteProfiles, ...)
  5  Signature accession (e.g. PF00069)
  6  Signature description
  7  Start location
  8  Stop location
  9  Score
  10 Status (T/? - almost always T)
  11 Date
  12 InterPro annotations - accession (e.g. IPR000719) or "-" if none
  13 InterPro annotations - description
  14 GO annotations (if -goterms was used)
  15 Pathway annotations (if -pa was used)

One protein can have many rows (one per domain hit, across analyses) - this collapses to
one row per protein, listing every distinct InterPro accession and every distinct Pfam
accession found for it (semicolon-separated), for
update_interpro_in_gff.py to consume.

Usage:
    python3 -m pangenome_pipeline.step04_interproscan.merge_interpro_results \\
        --results-dir interproscan_results/GCA_000002495.2 \\
        --out interpro_lookup/GCA_000002495.2_interpro.tsv
"""

import argparse
import glob
import os
import sys
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--results-dir", required=True, help="Directory with this genome's chunk_*_interproscan.tsv files"
    )
    parser.add_argument("--out", required=True, help="Output TSV: protein_id, interpro_ids, pfam_ids, n_domain_hits")
    parser.add_argument(
        "--fasta",
        default=None,
        help="Optional: this genome's own .faa - if given, reports every protein in it "
        "that never shows up in the combined TSV results. InterProScan's '-f tsv' "
        "format only ever writes a row when an analysis finds a hit, so a protein "
        "absent from every chunk's TSV most likely just has zero domain hits across "
        "all analyses run - completely normal, not a failure. This flag can't tell "
        "that apart from a genuinely malformed/truncated sequence InterProScan choked "
        "on, so treat its output as a list worth a manual spot-check, not a list of "
        "confirmed errors.",
    )
    args = parser.parse_args()

    tsv_files = sorted(glob.glob(os.path.join(args.results_dir, "chunk_*_interproscan.tsv")))
    if not tsv_files:
        sys.exit(f"ERROR: no chunk_*_interproscan.tsv files found in {args.results_dir}")

    # Sanity check: warn if a chunk_N.fasta-implied chunk has no matching TSV (InterProScan
    # failed/never ran on it) - compares chunk numbers seen in TSV filenames against the
    # full expected range, so a silently-missing chunk doesn't just produce an
    # incomplete-but-unflagged lookup table.
    chunk_nums_seen = set()
    for f in tsv_files:
        base = os.path.basename(f)
        try:
            n = int(base.split("_")[1])
            chunk_nums_seen.add(n)
        except (IndexError, ValueError):
            pass
    if chunk_nums_seen:
        expected = set(range(1, max(chunk_nums_seen) + 1))
        missing = expected - chunk_nums_seen
        if missing:
            print(
                f"  WARNING: chunk(s) {sorted(missing)} have no TSV in {args.results_dir} "
                f"- InterProScan may have failed on them. Lookup table will be incomplete "
                f"for those proteins.",
                file=sys.stderr,
            )

    interpro_ids = defaultdict(set)
    pfam_ids = defaultdict(set)
    n_hits = defaultdict(int)
    all_proteins = set()

    for tsv_path in tsv_files:
        with open(tsv_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 12:
                    continue
                protein_id = fields[0]
                analysis = fields[3]
                signature_acc = fields[4]
                interpro_acc = fields[11]

                all_proteins.add(protein_id)
                n_hits[protein_id] += 1

                if interpro_acc and interpro_acc != "-":
                    interpro_ids[protein_id].add(interpro_acc)
                if analysis.lower() == "pfam" and signature_acc and signature_acc != "-":
                    pfam_ids[protein_id].add(signature_acc)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out_fh:
        out_fh.write("protein_id\tinterpro_ids\tpfam_ids\tn_domain_hits\n")
        for protein_id in sorted(all_proteins):
            ipr = ";".join(sorted(interpro_ids[protein_id]))
            pfam = ";".join(sorted(pfam_ids[protein_id]))
            out_fh.write(f"{protein_id}\t{ipr}\t{pfam}\t{n_hits[protein_id]}\n")

    n_with_ipr = sum(1 for p in all_proteins if interpro_ids[p])
    n_with_pfam = sum(1 for p in all_proteins if pfam_ids[p])
    print(
        f"[{args.results_dir}] {len(all_proteins)} proteins scanned ({len(tsv_files)} chunk files) -> {args.out}",
        file=sys.stderr,
    )
    print(
        f"  {n_with_ipr} with >=1 InterPro ID ({100 * n_with_ipr / len(all_proteins):.1f}%), "
        f"{n_with_pfam} with >=1 Pfam ID ({100 * n_with_pfam / len(all_proteins):.1f}%)",
        file=sys.stderr,
    )

    if args.fasta:
        fasta_proteins = set()
        with open(args.fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    fasta_proteins.add(line[1:].split()[0])
        missing = fasta_proteins - all_proteins
        extra = all_proteins - fasta_proteins
        if not missing and not extra:
            print(
                f"  Completeness check: OK - all {len(fasta_proteins)} proteins in "
                f"{args.fasta} are present in the combined results.",
                file=sys.stderr,
            )
        else:
            if missing:
                print(
                    f"  NOTE: {len(missing)} / {len(fasta_proteins)} proteins in {args.fasta} "
                    f"never show up in any chunk TSV - almost certainly just proteins with zero "
                    f"domain hits across every analysis run (InterProScan's TSV format only writes "
                    f"a row per hit, none if there's no hit at all), not a processing failure. Not "
                    f"the same as a missing chunk *file*, which the check above already covers. "
                    f"First few, if you want to spot-check one: {sorted(missing)[:10]}",
                    file=sys.stderr,
                )
            if extra:
                print(
                    f"  NOTE: {len(extra)} protein IDs in the TSV results aren't in {args.fasta} - "
                    f"worth double-checking this is the right FASTA for this genome.",
                    file=sys.stderr,
                )


if __name__ == "__main__":
    main()
