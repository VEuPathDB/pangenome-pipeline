#!/usr/bin/env python3
"""
Check the ORF of every protein-coding gene of the final GFFs against its genome sequence, the same way for every gene
(Liftoff only checks the genes it transferred; a genome's own genes have no ORF result).

For each transcript, the CDS segments are joined in genome order (reverse-complemented on the minus strand) and checked:
  start_ok        starts with ATG
  stop_ok         ends with a stop codon (TAA, TAG, TGA)
  frame_ok        length is a multiple of 3
  internal_stops  stop codons before the last codon (standard genetic code)
A transcript has a valid ORF when start_ok, stop_ok and frame_ok are true and there is no internal stop. A gene has a
valid ORF when at least one of its transcripts has one. Genes with a pseudogene biotype are reported as such (their
ORF is still checked, for information).

Output: one TSV, one row per gene: genome, gene_id, original_ID, source, biotype, n_transcripts, best transcript,
cds_bp, start_ok, stop_ok, frame_ok, internal_stops, valid_orf.

With --out-gff-dir, every GFF is also written again with the result, in Liftoff's format: valid_ORFs=<number of
transcripts with a valid ORF> on the gene line, valid_ORF=True/False on each transcript line. An existing value is kept
when it agrees, replaced when it differs, and added when missing (the genome's own genes). The input GFFs are not changed.

Usage:
    python3 check_orf.py --gff-dir filtered_annotation_domain_only --gff-suffix _final_annotation_filtered.gff \
        --genome-dir genomes --genome-suffix .fna --genomes genome_order.txt -o orf_check/orf_check.tsv

Needs only the standard library.
"""
import argparse
import os
import re
import sys
from collections import defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
ORIG_RE = re.compile(r"(?:^|;)original_ID=([^;]+)")
BIOTYPE_RE = re.compile(r"(?:^|;)(?:ebi_biotype|gene_biotype|biotype)=([^;]+)")
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    seqs, name, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    return seqs


def read_gff(path):
    """gene -> {contig, strand, source, original_ID, biotype, transcripts: {tx: [(start, end)]}}"""
    parent, top, cds = {}, {}, defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
            if i:
                parent[i.group(1)] = p.group(1).split(",")[0] if p else None
                if not p:
                    o, b = ORIG_RE.search(f[8]), BIOTYPE_RE.search(f[8])
                    top[i.group(1)] = {"contig": f[0], "strand": f[6], "source": f[1],
                                       "original_ID": o.group(1) if o else i.group(1),
                                       "biotype": b.group(1) if b else ""}
            if f[2] == "CDS" and p:
                cds[p.group(1).split(",")[0]].append((int(f[3]), int(f[4])))

    def root(x):
        seen = set()
        while parent.get(x) and x not in seen:
            seen.add(x)
            x = parent[x]
        return x

    genes = {}
    for tx, segs in cds.items():
        g = root(tx)
        if g in top:
            genes.setdefault(g, dict(top[g], transcripts={}))["transcripts"][tx] = sorted(set(segs))
    return genes


def check(seq):
    n = len(seq) - len(seq) % 3
    codons = [seq[i:i + 3] for i in range(0, n, 3)]
    start_ok = seq[:3] == "ATG"
    stop_ok = len(seq) >= 3 and seq[-3:] in STOPS
    frame_ok = len(seq) % 3 == 0
    internal = sum(1 for c in codons[:-1] if c in STOPS)
    # a CDS with an N (assembly gap) is not a confirmed ORF, even when the gap keeps the frame
    no_gap = "N" not in seq.upper()
    return start_ok, stop_ok, frame_ok, internal, (start_ok and stop_ok and frame_ok and internal == 0 and no_gap)


def set_attr(field9, key, value):
    """Set key=value in a GFF attribute column; returns (new column, 'kept'|'changed'|'added')."""
    parts = field9.split(";")
    for i, kv in enumerate(parts):
        if kv.startswith(key + "="):
            if kv == f"{key}={value}":
                return field9, "kept"
            parts[i] = f"{key}={value}"
            return ";".join(parts), "changed"
    return field9.rstrip(";") + f";{key}={value}", "added"


def write_gff(src, dst, gene_valid, tx_valid):
    """Copy src to dst with valid_ORFs on gene lines and valid_ORF on transcript lines; returns counts per action."""
    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    counts = {"kept": 0, "changed": 0, "added": 0}
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                i = ID_RE.search(f[8])
                fid = i.group(1) if i else None
                if fid in gene_valid and "Parent=" not in f[8]:
                    f[8], act = set_attr(f[8], "valid_ORFs", gene_valid[fid])
                    counts[act] += 1
                elif fid in tx_valid:
                    f[8], act = set_attr(f[8], "valid_ORF", "True" if tx_valid[fid] else "False")
                    counts[act] += 1
                fout.write("\t".join(f) + "\n")
            else:
                fout.write(line)
    return counts["kept"], counts["changed"], counts["added"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--genomes", required=True, help="file with one genome per line, or comma-separated names")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--out-gff-dir", default="", help="also write every GFF with valid_ORFs / valid_ORF set from this check")
    args = ap.parse_args()

    if os.path.exists(args.genomes):
        genomes = [l.split()[0] for l in open(args.genomes) if l.strip() and not l.startswith("#")]
    else:
        genomes = [g.strip() for g in args.genomes.split(",") if g.strip()]

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        out.write("genome\tgene_id\toriginal_ID\tsource\tbiotype\tn_transcripts\tbest_transcript\tcds_bp\t"
                  "start_ok\tstop_ok\tframe_ok\tinternal_stops\tvalid_orf\n")
        for g in genomes:
            gff = os.path.join(args.gff_dir, g + args.gff_suffix)
            fa = os.path.join(args.genome_dir, g + args.genome_suffix)
            for p in (gff, fa):
                if not os.path.exists(p):
                    sys.exit(f"ERROR: missing {p}")
            seqs = read_fasta(fa)
            genes = read_gff(gff)
            n_valid = n_missing = 0
            tx_valid = {}
            for gid, rec in genes.items():
                contig = seqs.get(rec["contig"])
                if contig is None:
                    n_missing += 1
                    continue
                best = None
                for tx, segs in rec["transcripts"].items():
                    s = "".join(contig[a - 1:b] for a, b in segs)
                    if rec["strand"] == "-":
                        s = s.translate(COMP)[::-1]
                    r = check(s)
                    tx_valid[tx] = r[4]
                    key = (r[4], -r[3], r[0] + r[1] + r[2])
                    if best is None or key > best[0]:
                        best = (key, tx, len(s), r)
                _, tx, bp, (st, sp, fr, internal, valid) = best
                n_valid += valid
                out.write(f"{g}\t{gid}\t{rec['original_ID']}\t{rec['source']}\t{rec['biotype']}\t{len(rec['transcripts'])}\t"
                          f"{tx}\t{bp}\t{st}\t{sp}\t{fr}\t{internal}\t{valid}\n")
            msg = (f"{g}: {len(genes)} coding genes, {n_valid} with a valid ORF ({100 * n_valid / max(len(genes), 1):.1f}%)"
                   + (f"; {n_missing} on contigs missing from the FASTA" if n_missing else ""))
            if args.out_gff_dir:
                gene_valid = {gid: sum(1 for tx in rec["transcripts"] if tx_valid.get(tx)) for gid, rec in genes.items()
                              if any(tx in tx_valid for tx in rec["transcripts"])}
                kept, changed, added = write_gff(gff, os.path.join(args.out_gff_dir, g + args.gff_suffix), gene_valid, tx_valid)
                msg += f"; GFF attributes kept {kept}, changed {changed}, added {added}"
            print(msg, file=sys.stderr)
    print(f"Wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
