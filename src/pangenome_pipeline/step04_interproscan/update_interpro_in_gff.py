#!/usr/bin/env python3
"""
Writes this genome's own step-04 InterProScan results into its own GFF, replacing whatever
InterPro/PFAM evidence (if any) the original annotation already carried - so every genome ends
up on the same InterProScan version, whether its original submitter ran one, ran an older one,
or never ran one at all.

This is deliberately simple, one genome at a time, no cross-genome logic: at this point in the
pipeline every gene's protein still belongs to its own genome (Liftoff/all-x-all hasn't run
yet), so there's no "which genome did this gene actually come from" question to resolve - that
version of the problem belongs to step 05's own injection, later, after transfer.

Matching: AGAT (`agat_sp_extract_sequences.pl -p -t cds`, run in extract_proteins.slurm.sh)
uses each mRNA feature's own ID= as the header of its translated protein - that's the
`protein_id` merge_interpro_results.py's lookup table is keyed by. So for each CDS line, its
own Parent= (the mRNA it belongs to) is used as the lookup key directly - no separate ID
resolution needed.

For each CDS line's Dbxref= attribute:
  1. Strip any existing "InterPro:..."/"PFAM:..." entries (old evidence, any source).
  2. Look up this CDS's mRNA (its Parent=) in the lookup table.
  3. If found, append fresh "InterPro:<id>"/"PFAM:<id>" entries, one token per ID.
Everything else already in Dbxref= (e.g. "NCBI_GP:...") is left untouched. A CDS whose Dbxref=
becomes empty after stripping, with nothing new to add, has the attribute removed entirely.

Usage:
    python3 -m pangenome_pipeline.step04_interproscan.update_interpro_in_gff --gff gff_cleaned/AF293.gff3 \\
        --lookup interpro/interpro_lookup/AF293_interpro.tsv \\
        --out gff_with_interpro/AF293.gff3
"""

import argparse
import re
import sys

from pangenome_pipeline.common.gff import DBXREF_RE, PARENT_RE


def load_lookup(path):
    table = {}
    with open(path) as fh:
        fh.readline()  # skip the header row
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            protein_id, interpro_ids, pfam_ids = fields[0], fields[1], fields[2]
            table[protein_id] = (interpro_ids, pfam_ids)
    return table


def strip_old_domain_xrefs(dbxref_value):
    """Remove InterPro:/PFAM: entries (any case) from a Dbxref= value, keep the rest."""
    kept = [tok for tok in dbxref_value.split(",") if not tok.strip().lower().startswith(("interpro:", "pfam:"))]
    return kept


def new_domain_tokens(interpro_ids, pfam_ids):
    tokens = []
    if interpro_ids:
        tokens += [f"InterPro:{i}" for i in interpro_ids.split(";") if i]
    if pfam_ids:
        tokens += [f"PFAM:{i}" for i in pfam_ids.split(";") if i]
    return tokens


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--lookup", required=True, help="this genome's own <genome>_interpro.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lookup = load_lookup(args.lookup)

    with open(args.gff) as fh:
        raw_lines = fh.readlines()

    n_cds = 0
    n_had_old = 0
    n_matched = 0
    n_no_mrna_id = 0

    out_lines = []
    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            out_lines.append(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "CDS":
            out_lines.append(line)
            continue

        n_cds += 1
        attrs = fields[8]
        p = PARENT_RE.search(attrs)
        mrna_id = p.group(1) if p else None
        if not mrna_id:
            n_no_mrna_id += 1
            out_lines.append(line)
            continue

        m = DBXREF_RE.search(attrs)
        kept_tokens = strip_old_domain_xrefs(m.group(1)) if m else []
        if m and len(kept_tokens) < len(m.group(1).split(",")):
            n_had_old += 1

        hit = lookup.get(mrna_id)
        if hit:
            interpro_ids, pfam_ids = hit
            added = new_domain_tokens(interpro_ids, pfam_ids)
            if added:
                n_matched += 1
                kept_tokens += added

        new_dbxref = ",".join(kept_tokens)
        if m:
            old_attr = f"Dbxref={m.group(1)}"
            if new_dbxref:
                attrs = attrs.replace(old_attr, f"Dbxref={new_dbxref}", 1)
            else:
                # drop the attribute entirely, plus one leading/trailing ';' if left dangling
                attrs = re.sub(r";?Dbxref=[^;]+", "", attrs, count=1)
        elif new_dbxref:
            attrs = attrs + f";Dbxref={new_dbxref}"
        fields[8] = attrs
        out_lines.append("\t".join(fields) + "\n")

    with open(args.out, "w") as out_fh:
        out_fh.writelines(out_lines)

    pct = 100 * n_matched / n_cds if n_cds else 0
    print(
        f"[{args.gff}] {n_cds} CDS line(s): {n_had_old} had old InterPro/PFAM evidence "
        f"stripped, {n_matched} ({pct:.1f}%) got fresh evidence -> {args.out}",
        file=sys.stderr,
    )
    if n_no_mrna_id:
        print(f"  NOTE: {n_no_mrna_id} CDS line(s) had no Parent= - left untouched.", file=sys.stderr)


if __name__ == "__main__":
    main()
