#!/usr/bin/env python3
"""
Trim FASTA headers and normalise the filename. No other input needed.

For every FASTA given, it:

  1. cuts each header down to its first token
         >CM000171.1 Aspergillus fumigatus Af293 chromosome 1   ->   >CM000171.1
     Nothing is removed - every contig is kept, only the text after the first
     token is dropped. Those trimmed IDs are what must match column 1 of the GFF.

  2. writes it out as <basename><ext>, where <basename> is
       - the NCBI accession, when the name starts with one:
             GCA_040142875.1_UCR_Afum_08-19-02-30_1.0_genomic.fna  ->  GCA_040142875.1.fna
       - otherwise the name unchanged.

Usage
-----
# in place (typical)
rename_headers.py genomes/*.fna genomes/*.fasta --in-place

# keep the originals - write trimmed copies into another folder, all as .fna
rename_headers.py raw/*.fna --outdir genomes --to-ext .fna

# just check the trimmed IDs against the GFFs, touch nothing
rename_headers.py genomes/*.fna --check-gff gff_originals
"""

import argparse
import os
import re
import sys
from functools import partial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.gff import find_gff as _find_gff
from pangenome_common.gff import gff_seqids

# Also accepts gzipped GFFs, unlike most steps.
find_gff = partial(_find_gff, exts=(".gff3", ".gff", ".gff3.gz", ".gff.gz"))

FASTA_EXTS = (".fna", ".fasta", ".fa")
ACC_RE = re.compile(r"^(GC[AF]_\d+\.\d+)")


def stem_and_ext(path):
    base = os.path.basename(path)
    for ext in FASTA_EXTS:
        if base.lower().endswith(ext):
            return base[: -len(ext)], base[-len(ext):]
    root, ext = os.path.splitext(base)
    return root, ext


def basename_for(stem):
    """NCBI accession if the name starts with one, else the stem unchanged."""
    m = ACC_RE.match(stem)
    return m.group(1) if m else stem


def trim_headers(in_path, out_path):
    """Trim headers to their first token. out_path=None -> read the IDs, write nothing."""
    n_seqs = 0
    seen = set()
    dupes = set()
    ids = []
    outfile = open(out_path, "w") if out_path is not None else None
    try:
        with open(in_path) as infile:
            for line in infile:
                if line.startswith(">"):
                    n_seqs += 1
                    header = line[1:].strip().split()[0] if line[1:].strip() else ""
                    if not header:
                        sys.exit(f"ERROR: empty FASTA header in {in_path}")
                    if header in seen:
                        dupes.add(header)
                    seen.add(header)
                    ids.append(header)
                    if outfile:
                        outfile.write(f">{header}\n")
                elif outfile:
                    outfile.write(line)
    finally:
        if outfile:
            outfile.close()
    if dupes:
        if out_path is not None and out_path != in_path:
            try:
                os.remove(out_path)
            except OSError:
                pass
        sys.exit(f"ERROR: {in_path}: headers collide after trimming to first token: "
                 f"{sorted(dupes)[:10]} - fix the source assembly.")
    return n_seqs, ids


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="FASTA file(s): .fna / .fasta / .fa")
    ap.add_argument("--in-place", action="store_true",
                    help="edit the files where they are")
    ap.add_argument("--outdir", metavar="DIR",
                    help="write trimmed copies here instead (keeps the originals)")
    ap.add_argument("--to-ext", choices=FASTA_EXTS,
                    help="force every output to this extension (default: keep the input's)")
    ap.add_argument("--check-gff", metavar="DIR",
                    help="don't edit - just report any sequence ID missing from "
                         "DIR/<basename>.gff3 column 1")
    args = ap.parse_args()

    check_only = bool(args.check_gff) and not args.in_place and not args.outdir
    if args.in_place and args.outdir:
        sys.exit("ERROR: use --in-place or --outdir, not both.")
    if not (args.in_place or args.outdir or check_only):
        sys.exit("ERROR: give --in-place, --outdir, or --check-gff.")

    # Resolve output names up front and catch collisions before touching anything.
    plan = []          # (in_path, stem, ext, out_ext, base)
    errors = []
    for in_path in args.inputs:
        if not os.path.isfile(in_path):
            errors.append(f"no such file: {in_path}")
            continue
        stem, ext = stem_and_ext(in_path)
        plan.append((in_path, stem, ext, args.to_ext or ext, basename_for(stem)))

    seen = {}
    for in_path, _, _, out_ext, base in plan:
        key = base + out_ext
        if key in seen:
            errors.append(f"{in_path} and {seen[key]} would both become {key}")
        else:
            seen[key] = in_path
    if errors:
        sys.exit("ERROR:\n  " + "\n  ".join(errors))

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    problems = 0
    for in_path, stem, ext, out_ext, base in plan:
        if check_only:
            n_seqs, ids = trim_headers(in_path, None)
            print(f"{in_path}   ({n_seqs} sequences)")
        elif args.in_place:
            tmp = in_path + ".tmp_rename"
            n_seqs, ids = trim_headers(in_path, tmp)
            os.replace(tmp, in_path)
            out_path = in_path
            if base + out_ext != os.path.basename(in_path):
                out_path = os.path.join(os.path.dirname(in_path), base + out_ext)
                os.replace(in_path, out_path)
            print(f"{in_path}  ->  {out_path}   ({n_seqs} sequences)")
        else:
            out_path = os.path.join(args.outdir, base + out_ext)
            n_seqs, ids = trim_headers(in_path, out_path)
            print(f"{in_path}  ->  {out_path}   ({n_seqs} sequences)")

        if args.check_gff:
            gff = find_gff(args.check_gff, base) or find_gff(args.check_gff, stem)
            if not gff:
                print(f"  check-gff: no GFF for '{base}' in {args.check_gff} - skipped",
                      file=sys.stderr)
            else:
                missing = [i for i in ids if i not in gff_seqids(gff)]
                if missing:
                    problems += 1
                    print(f"  WARNING: {len(missing)}/{len(ids)} trimmed IDs are not in "
                          f"{os.path.basename(gff)} col 1 (first: {missing[:5]})",
                          file=sys.stderr)
                else:
                    print(f"  check-gff: OK - all {len(ids)} IDs in {os.path.basename(gff)}")

    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
