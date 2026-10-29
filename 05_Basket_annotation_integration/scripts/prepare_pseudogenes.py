#!/usr/bin/env python3
"""
First step of the basket annotation: makes the pseudogenes usable by the next steps, which work on CDS.

A pseudogene of the annotations (a gene with a pseudogenic_transcript, or ebi_biotype=pseudogene, and no CDS) gets:
  - its transcript typed "transcript" (pseudogenic_transcript is not read by the next steps nor by get_pangenes);
  - one CDS line per transcript, from the transcript's first to its last base (its span, introns and gaps included:
    the real sequence of the locus, defect included), tagged pseudogene_span=True; overlap_filter.py removes it
    again from the final GFFs, so a pseudogene is written as gene + transcript + exon, as in the source;
  - pseudogene=True on its gene line (typed protein_coding_gene when it was typed pseudogene).
Only a gene annotated as a pseudogene in its own genome is called a pseudogene. With --transfers (the Liftoff
transfers), a pseudogene's copy in another genome gets the span CDS too (tagged cds_from_transcript_span=True, kept
in the final GFFs) but not pseudogene=True: it is a transferred model like any other (its ORF checked, miniprot).
Only the reference's pseudogenes are used (--reference): they win their basket like every reference gene. Those of
the other genomes are left as they are (without a CDS, the next steps do not use them): a singleton pseudogene adds
nothing to the pangenome, and where another genome has the gene, its transfer covers the place.
Every other line is copied unchanged. A GFF without such pseudogenes is linked, not copied (the transfer folder is
large).

Usage: python3 prepare_pseudogenes.py --in-dir gff_for_liftoff --out-dir basket_annotation/input/native
Writes the GFFs to --out-dir and prints the number of pseudogenes per file.
"""
import argparse
import glob
import os
import re
from collections import defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def convert(path, own):
    """the lines of path with its pseudogenes made usable; None when it has none."""
    lines = open(path).read().splitlines()
    rows = [l.split("\t") if not l.startswith("#") else None for l in lines]
    gene_attr, tx_of, tx_type, has_cds = {}, {}, {}, set()
    for f in rows:
        if not f or len(f) < 9:
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if i and not p:
            gene_attr[i.group(1)] = f[8]
        elif i and p and f[2] in ("mRNA", "transcript", "pseudogenic_transcript"):
            tx_of[i.group(1)] = p.group(1).split(",")[0]
            tx_type[i.group(1)] = f[2]
        elif f[2] == "CDS" and p:
            for par in p.group(1).split(","):
                has_cds.add(tx_of.get(par, par))
    with_ptx = {gg for t, gg in tx_of.items() if tx_type[t] == "pseudogenic_transcript"}
    pseudo = {g for g, attr in gene_attr.items()
              if g not in has_cds and (g in with_ptx or re.search(r"(?:^|;)ebi_biotype=pseudogene(?:;|$)", attr))}
    if not pseudo:
        return None, 0
    out = []
    for line, f in zip(lines, rows):
        if not f or len(f) < 9:
            out.append(line)
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if i and not p and i.group(1) in pseudo:
            typ = "protein_coding_gene" if f[2] == "pseudogene" else f[2]
            out.append("\t".join(f[:2] + [typ] + f[3:8] + [f[8].rstrip(";") + (";pseudogene=True" if own else "")]))
        elif i and p and i.group(1) in tx_of and tx_of[i.group(1)] in pseudo:
            tid = i.group(1)
            out.append("\t".join(f[:2] + ["transcript"] + f[3:]))
            tag = "pseudogene_span=True" if own else "cds_from_transcript_span=True"
            out.append("\t".join([f[0], f[1], "CDS", f[3], f[4], ".", f[6], "0", f"ID={tid}-span_cds;Parent={tid};{tag}"]))
        else:
            out.append(line)
    return out, len(pseudo)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--pattern", default="*", help="files to take from --in-dir (default: all GFFs)")
    ap.add_argument("--reference", required=True,
                    help="the reference genome: only its pseudogenes are used (its own GFF, and its transfers "
                         "<reference>_onto_*); the other genomes' pseudogenes are left out, as they are")
    ap.add_argument("--transfers", action="store_true",
                    help="--in-dir holds Liftoff transfers: a pseudogene's copy in another genome is not called a "
                         "pseudogene (no pseudogene=True) and keeps its span as CDS, like any transferred model")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    total = defaultdict(int)
    for path in sorted(glob.glob(os.path.join(a.in_dir, a.pattern))):
        name = os.path.basename(path)
        if not re.search(r"\.gff3?$|_polished$", name):
            continue
        dst = os.path.join(a.out_dir, name)
        if os.path.lexists(dst):
            os.remove(dst)
        mine = name.startswith(f"{a.reference}_onto_") if a.transfers else re.sub(r"\.gff3?$", "", name) == a.reference
        out, n = convert(path, own=not a.transfers) if mine else (None, 0)
        if out is None:
            os.symlink(os.path.realpath(path), dst)
        else:
            with open(dst, "w") as fh:
                fh.write("\n".join(out) + "\n")
            print(f"  {name}: {n} pseudogenes")
        total["files"] += 1
        total["pseudogenes"] += n
    print(f"{total['files']} GFFs; {total['pseudogenes']} pseudogenes made usable")


if __name__ == "__main__":
    main()
