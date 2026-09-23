#!/usr/bin/env python3
"""
Extract IDs of "good transfer" protein_coding_gene features from a Liftoff GFF3.

A gene is kept when ALL of the following hold:
  - it has a valid ORF (valid_ORF=True, or valid_ORFs >= 1)
  - alignment coverage >= --min-coverage (default 0.9)
  - sequence identity   >= --min-identity (default 0.9)
  - it is the primary copy, i.e. extra_copy_number == 0 (or absent)
    (use --keep-extra-copies to keep extra copies too)

Output is a plain text file with one gene ID per line, suitable for
agat_sp_filter_feature_from_keep_list.pl's --keep_list.

Usage:
    python3 get_good_transfer_ids.py --gff STRAIN_liftoff_updated.gff3_polished --out keep_ids.txt
"""

import argparse
import re
import sys

# Anchored on the start of the column or a ";" so a name like "ref_sequence_ID=" or "xID=" is never
# mistaken for the attribute we want. Numbers are read as whole values (anything up to the next ";")
# and converted with float(), so e.g. "1e-05" is not silently read as "1".
ATTR_PATTERNS = {
    "ID": re.compile(r"(?:^|;)ID=([^;]+)"),
    "valid_ORF": re.compile(r"(?:^|;)valid_ORF=(True|False)"),
    "valid_ORFs": re.compile(r"(?:^|;)valid_ORFs=(\d+)"),
    "coverage": re.compile(r"(?:^|;)coverage=([^;]+)"),
    "sequence_ID": re.compile(r"(?:^|;)sequence_ID=([^;]+)"),
    "extra_copy_number": re.compile(r"(?:^|;)extra_copy_number=(\d+)"),
}


def to_float(value):
    """A missing or unparsable number counts as 0 - i.e. never 'good'."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def parse_attributes(attribute_field):
    """Pull out the fields we need from a GFF attribute column."""
    parsed = {}
    for key, pattern in ATTR_PATTERNS.items():
        match = pattern.search(attribute_field)
        parsed[key] = match.group(1) if match else None
    return parsed


def has_valid_orf(attrs):
    """Prefer the explicit valid_ORF=True/False flag (mRNA-level); fall back
    to valid_ORFs=N (gene-level count) if that's what's present instead."""
    if attrs["valid_ORF"] == "True":
        return True
    if attrs["valid_ORF"] == "False":
        return False
    if attrs["valid_ORFs"] is not None:
        return int(attrs["valid_ORFs"]) >= 1
    return False


def is_good_transfer(attrs, min_coverage, min_identity, keep_extra_copies):
    if not has_valid_orf(attrs):
        return False

    coverage = to_float(attrs["coverage"])
    identity = to_float(attrs["sequence_ID"])
    if coverage < min_coverage or identity < min_identity:
        return False

    if not keep_extra_copies:
        extra_copy_number = int(attrs["extra_copy_number"]) if attrs["extra_copy_number"] is not None else 0
        if extra_copy_number > 0:
            return False

    return True


def extract_good_transfer_ids(gff_path, min_coverage, min_identity, keep_extra_copies, gene_feature_types):
    ids = []
    n_protein_coding_genes = 0
    with open(gff_path) as gff_file:
        for line in gff_file:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in gene_feature_types:
                continue
            n_protein_coding_genes += 1
            attrs = parse_attributes(fields[8])
            if attrs["ID"] is None:
                continue
            if is_good_transfer(attrs, min_coverage, min_identity, keep_extra_copies):
                ids.append(attrs["ID"])
    return ids, n_protein_coding_genes


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gff", required=True, help="Input Liftoff GFF3 (the *_strain_liftoff_updated.gff3_polished file)")
    parser.add_argument("--out", required=True, help="Output text file, one gene ID per line")
    parser.add_argument("--min-coverage", type=float, default=0.9, help="Minimum Liftoff coverage (default: 0.9)")
    parser.add_argument("--min-identity", type=float, default=0.9, help="Minimum Liftoff sequence identity (default: 0.9)")
    parser.add_argument("--keep-extra-copies", action="store_true",
                         help="Keep extra copies (extra_copy_number > 0) too, instead of restricting to the primary copy")
    parser.add_argument("--gene-feature-types", default="gene,protein_coding_gene",
                         help="Comma-separated level-1 feature type name(s) that mark a protein-coding gene in "
                              "the input GFF (default: 'gene,protein_coding_gene', covering both common "
                              "conventions - NCBI-style GFFs typically use 'gene', Liftoff/VEuPathGDB-style "
                              "references typically use 'protein_coding_gene')")
    args = parser.parse_args()

    gene_feature_types = {t.strip() for t in args.gene_feature_types.split(",") if t.strip()}

    ids, n_total = extract_good_transfer_ids(
        args.gff, args.min_coverage, args.min_identity, args.keep_extra_copies, gene_feature_types
    )

    with open(args.out, "w") as out_file:
        out_file.write("\n".join(ids) + "\n" if ids else "")

    print(f"[{args.gff}] {len(ids)} / {n_total} gene features kept -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
