#!/usr/bin/env python3
"""
Remove named contigs from the assemblies (and the matching features from the
GFFs). Use it to strip organelle contigs found by the organelle screen (section
7.2 / 1.2), or any other sequence you have decided not to carry into the
pangenome.

It is driven by a plain list file, one  <genome_basename>  <contig_id>  pair per
line (whitespace- or tab-separated); '#' comments and blank lines are ignored:

    # genome            contig to drop
    AF293               mito_A_fumigatus_Af293
    GCA_040126065.1     CP146368.1
    GCA_040142845.1     JBBFMX010000010.1
    GCA_040142845.1     JBBFMX010000011.1

For every genome named in the list it:

  1. rewrites  <genomes-dir>/<genome>.fna  without those sequences, and
  2. rewrites  <gff-dir>/<genome>.gff3     without the lines whose column 1 is a
     dropped contig  (also  ##sequence-region  pragmas and any  ##FASTA  records).

Contig IDs are matched against the FASTA header's first token (the same trimmed
ID that must match GFF column 1 - see section 5.1). A contig already gone from the
FASTA but still in the GFF is fine - it just cleans the GFF (so you can run this
after stripping the assemblies some other way). A contig in neither the FASTA nor
the GFF is skipped with a note - so a cumulative drop list stays safe to re-run
after some of its contigs have already been removed. Pass  --strict  to make that
an error instead (catches typos).

Usage
-----
# dry run - the default: report what would be removed, touch nothing
drop_contigs.py drop_contigs.txt

# do it: rewrite genomes/ and gff_originals/ in place
drop_contigs.py drop_contigs.txt --apply

# assemblies only, no GFFs
drop_contigs.py drop_contigs.txt --apply --gff-dir ""
"""

import argparse
import os
import sys
from collections import defaultdict

FASTA_EXTS = (".fna", ".fasta", ".fa")
GFF_EXTS = (".gff3", ".gff")


def read_list(path):
    """-> {genome: set(contig_ids)}, preserving nothing else."""
    want = defaultdict(set)
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                sys.exit(f"ERROR: {path}:{n}: expected '<genome> <contig>', got: {line!r}")
            genome, contig = parts
            want[genome].add(contig)
    if not want:
        sys.exit(f"ERROR: no (genome, contig) pairs in {path}")
    return want


def find_one(directory, base, exts):
    for e in exts:
        cand = os.path.join(directory, base + e)
        if os.path.isfile(cand):
            return cand
    return None


def fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                tok = line[1:].strip().split()
                ids.append(tok[0] if tok else "")
    return ids


def gff_seqids(path):
    """Every sequence name the GFF refers to: column 1, ##sequence-region, ##FASTA headers."""
    ids = set()
    in_fasta = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if in_fasta:
                if line.startswith(">"):
                    tok = line[1:].strip().split()
                    if tok:
                        ids.add(tok[0])
                continue
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 2:
                    ids.add(parts[1])
                continue
            if line.startswith("#") or not line.strip():
                continue
            ids.add(line.split("\t", 1)[0])
    return ids


def filter_fasta(in_path, drop, out_fh):
    """Write in_path to out_fh without records whose first-token ID is in `drop`.
    Returns {dropped_id: length_bp}."""
    dropped = {}
    keep = True
    cur_id = None
    cur_len = 0
    with open(in_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id is not None and not keep:
                    dropped[cur_id] = cur_len
                tok = line[1:].strip().split()
                cur_id = tok[0] if tok else ""
                cur_len = 0
                keep = cur_id not in drop
                if keep:
                    out_fh.write(line)
            else:
                if keep:
                    out_fh.write(line)
                else:
                    cur_len += len(line.strip())
        if cur_id is not None and not keep:
            dropped[cur_id] = cur_len
    return dropped


def filter_gff(in_path, drop, out_fh):
    """Copy in_path to out_fh, dropping annotation lines / pragmas / embedded FASTA
    records that belong to a contig in `drop`. Returns the number of feature lines removed."""
    removed = 0
    in_fasta = False
    keep_fasta_rec = True
    with open(in_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                out_fh.write(line)
                continue
            if in_fasta:
                if line.startswith(">"):
                    tok = line[1:].strip().split()
                    keep_fasta_rec = (tok[0] if tok else "") not in drop
                if keep_fasta_rec:
                    out_fh.write(line)
                continue
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 2 and parts[1] in drop:
                    continue
                out_fh.write(line)
                continue
            if line.startswith("#") or not line.strip():
                out_fh.write(line)
                continue
            if line.split("\t", 1)[0] in drop:
                removed += 1
                continue
            out_fh.write(line)
    return removed


def write_filtered(in_path, drop, filt):
    """Run `filt(in_path, drop, tmp_fh)` to a temp file, then atomically replace in_path."""
    tmp = in_path + ".tmp_dropcontigs"
    with open(tmp, "w") as out_fh:
        result = filt(in_path, drop, out_fh)
    os.replace(tmp, in_path)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("list", help="the <genome> <contig> list file")
    ap.add_argument(
        "--genomes-dir",
        default="genomes",
        metavar="DIR",
        help="assemblies as <genome>.fna / .fasta / .fa (default: genomes)",
    )
    ap.add_argument(
        "--gff-dir",
        default="gff_originals",
        metavar="DIR",
        help='annotations as <genome>.gff3 / .gff; "" to skip GFFs (default: gff_originals)',
    )
    ap.add_argument("--apply", action="store_true", help="actually rewrite the files (default: dry run - report only)")
    ap.add_argument(
        "--strict",
        action="store_true",
        help="fail if a listed contig is in neither the FASTA nor the GFF "
        "(default: skip it with a note - so a cumulative drop list stays "
        "safe to re-run after some contigs are already removed)",
    )
    args = ap.parse_args()

    want = read_list(args.list)

    # Resolve every file and validate every contig ID up front - touch nothing until
    # the whole plan is known to be sound.
    plan = []  # (genome, fna_path_or_None, gff_path_or_None, drop_set, in_fa_set)
    errors = []
    for genome in sorted(want):
        drop = want[genome]
        fna = find_one(args.genomes_dir, genome, FASTA_EXTS)
        gff = find_one(args.gff_dir, genome, GFF_EXTS) if args.gff_dir else None
        if not fna and not gff:
            errors.append(f"{genome}: no FASTA in {args.genomes_dir}/ and no GFF in {args.gff_dir or '(disabled)'}/")
            continue
        in_fa = set(fasta_ids(fna)) & drop if fna else set()
        in_gff = gff_seqids(gff) & drop if gff else set()
        unknown = sorted(drop - in_fa - in_gff)
        if unknown:
            where = os.path.basename(fna) if fna else ""
            where += (" / " if fna and gff else "") + (os.path.basename(gff) if gff else "")
            if args.strict:
                errors.append(f"{genome}: contig(s) not in {where}: {unknown}")
            else:
                # already removed by an earlier run (the list is a cumulative
                # manifest), or a typo - either way nothing to do. Skip, don't fail.
                print(
                    f"note: {genome}: {unknown} not in {where} (already removed, or a typo) - skipping", file=sys.stderr
                )
                drop = drop - set(unknown)
        only_gff = sorted(in_gff - in_fa)
        if only_gff and fna:
            print(
                f"note: {genome}: {only_gff} already absent from {os.path.basename(fna)} - will clean the GFF only",
                file=sys.stderr,
            )
        if not gff and args.gff_dir:
            print(f"note: {genome}: no GFF in {args.gff_dir}/ - assembly only", file=sys.stderr)
        if not drop:
            continue  # nothing left for this genome
        plan.append((genome, fna, gff, drop, in_fa))

    if errors:
        sys.exit("ERROR:\n  " + "\n  ".join(errors))
    if not plan:
        print("Nothing to do - every listed contig is already gone.")
        return

    verb = "Removing" if args.apply else "Would remove"
    for genome, fna, gff, drop, in_fa in plan:
        lens = {i: L for i, L in _fasta_lengths(fna).items() if i in drop} if fna else {}
        print(f"\n{genome}")
        for c in sorted(drop):
            tag = f"{lens[c]:,} bp" if c in in_fa else "not in FASTA - GFF only"
            print(f"  {verb}: {c}  ({tag})")
        if args.apply:
            if fna and in_fa:
                got = write_filtered(fna, drop, filter_fasta)
                kept = len(fasta_ids(fna))
                print(f"  {os.path.basename(fna)}: dropped {len(got)} contig(s), {kept} left")
            if gff:
                nrem = write_filtered(gff, drop, filter_gff)
                print(f"  {os.path.basename(gff)}: dropped {nrem} feature line(s)")
        elif gff:
            n = _gff_feature_count(gff, drop)
            print(f"  {os.path.basename(gff)}: {n} feature line(s) would go")

    if not args.apply:
        print("\n(dry run - nothing written; re-run with --apply)")


def _fasta_lengths(path):
    out, cur, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    out[cur] = n
                tok = line[1:].strip().split()
                cur, n = (tok[0] if tok else ""), 0
            else:
                n += len(line.strip())
    if cur is not None:
        out[cur] = n
    return out


def _gff_feature_count(path, drop):
    n = 0
    in_fasta = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
            if in_fasta or line.startswith("#") or not line.strip():
                continue
            if line.split("\t", 1)[0] in drop:
                n += 1
    return n


if __name__ == "__main__":
    main()
