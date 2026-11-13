#!/usr/bin/env python3
"""
Protein evidence for the genes the pipeline added to the reference (not_in_reference_annotation=True): is the protein
known in another species? Each protein is searched with NCBI blastp against nr over the internet (blastp -remote),
with the gene's own species left out (--exclude-taxid), so a hit is evidence from outside the species.

Writes into --out-dir (nothing else is touched):
  new_genes.faa            the proteins of the new genes (from --gff and --fasta)
  batches/batch_NNN.tsv    the raw blastp table of each batch of --batch proteins (outfmt 6, columns in BLAST_COLS)
  protein_evidence.tsv     one line per gene: its best hit outside the species, or "no hit"
A batch already in batches/ is not searched again, so a run that stops (time limit, network) is resumed by submitting
it again. --summary-only rebuilds protein_evidence.tsv from the batches done so far.

Usage:
    python3 protein_evidence.py --gff basket_annotation/gff_prefer_panel_model_final/AF293_basket_annotation.gff \\
        --fasta genomes/AF293.fna --exclude-taxid 746128 --out-dir reference_check/protein_evidence
"""
import argparse
import csv
import glob
import gzip
import os
import re
import subprocess
import sys
import time
from collections import defaultdict

att = lambda s, k: (re.search(r"(?:^|;)" + k + r"=([^;\n]*)", s) or [None, ""])[1]
openf = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {x + y + z: AA[16 * i + 4 * j + k] for i, x in enumerate(B) for j, y in enumerate(B) for k, z in enumerate(B)}
BLAST_COLS = "qseqid sacc pident length qlen slen evalue bitscore qcovs staxids sscinames stitle"


def new_genes(path):
    """new gene ID -> (contig, strand, [CDS (start, end)]) from the longest transcript"""
    genes, t2g, cds = {}, {}, defaultdict(list)
    for l in openf(path):
        if l.startswith("#") or l.count("\t") < 8:
            continue
        f = l.rstrip("\n").split("\t")
        if f[2] == "gene":
            if att(f[8], "not_in_reference_annotation") == "True":
                genes[att(f[8], "ID")] = (f[0], f[6])
        elif f[2] in ("mRNA", "transcript"):
            t2g[att(f[8], "ID")] = att(f[8], "Parent")
        elif f[2] == "CDS":
            for p in att(f[8], "Parent").split(","):
                cds[p].append((int(f[3]), int(f[4])))
    best = {}
    for t, segs in cds.items():
        g = t2g.get(t, t)
        if g in genes and len(segs) > len(best.get(g, [])):
            best[g] = sorted(segs)
    return {g: (genes[g][0], genes[g][1], best[g]) for g in genes if g in best}


def fasta(path, keep):
    seqs, name, parts = {}, None, []
    for l in openf(path):
        if l.startswith(">"):
            if name in keep:
                seqs[name] = "".join(parts)
            name, parts = l[1:].split()[0], []
        elif name in keep:
            parts.append(l.strip())
    if name in keep:
        seqs[name] = "".join(parts)
    return seqs


def protein(seq, segs, strand):
    s = "".join(seq[a - 1:b] for a, b in segs).upper()
    if strand == "-":
        s = s.translate(COMP)[::-1]
    p = "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
    p = p[:-1] if p.endswith("*") else p
    return p.replace("*", "X")   # a broken ORF: its internal stops as X, so blastp takes it


def species(row):
    """the species of a hit: sscinames when the remote search fills it, else the [Species] of the title"""
    if row.get("sscinames") and row["sscinames"] != "N/A":
        return row["sscinames"].split(";")[0]
    m = re.findall(r"\[([^\[\]]+)\]", row.get("stitle", ""))
    return m[-1] if m else ""


def summary(out_dir, ids):
    best = {}
    for p in sorted(glob.glob(os.path.join(out_dir, "batches", "batch_*.tsv"))):
        for l in open(p):
            f = l.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            r = dict(zip(BLAST_COLS.split(), f))
            r["bitscore"] = float(r["bitscore"])
            n = best.setdefault(r["qseqid"], [None, 0])
            n[1] += 1
            if n[0] is None or r["bitscore"] > n[0]["bitscore"]:
                n[0] = r
    done = set()
    for p in glob.glob(os.path.join(out_dir, "batches", "batch_*.ids")):
        done.update(open(p).read().split())
    out = os.path.join(out_dir, "protein_evidence.tsv")
    with open(out, "w") as w:
        w.write("gene\tresult\tbest_hit\tspecies\tidentity\tquery_coverage\tevalue\thits\ttitle\n")
        for g in ids:
            if g in best:
                r, n = best[g]
                w.write("\t".join([g, "hit", r["sacc"], species(r), r["pident"], r["qcovs"], r["evalue"], str(n),
                                   r["stitle"].replace("\t", " ")]) + "\n")
            elif g in done:
                w.write("\t".join([g, "no hit", "", "", "", "", "", "0", ""]) + "\n")
            else:
                w.write("\t".join([g, "not searched", "", "", "", "", "", "", ""]) + "\n")
    hit = sum(1 for g in ids if g in best)
    print(f"{len(ids)} genes: {hit} with a hit outside the species, {len(done & set(ids)) - hit} no hit, "
          f"{len(set(ids) - done)} not searched yet -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, help="the reference's final GFF (step 05)")
    ap.add_argument("--fasta", required=True, help="the reference genome (.fna, plain or .gz)")
    ap.add_argument("--exclude-taxid", required=True, help="the species left out of the search (e.g. 746128 A. fumigatus, 5811 T. gondii)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--db", default="nr")
    ap.add_argument("--evalue", default="1e-5")
    ap.add_argument("--batch", type=int, default=10, help="proteins per remote search")
    ap.add_argument("--pause", type=int, default=15, help="seconds between two remote searches (NCBI asks for >= 10)")
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args()

    os.makedirs(os.path.join(a.out_dir, "batches"), exist_ok=True)
    genes = new_genes(a.gff)
    ids = sorted(genes)
    faa = os.path.join(a.out_dir, "new_genes.faa")
    if not os.path.exists(faa):
        seqs = fasta(a.fasta, {c for c, _, _ in genes.values()})
        with open(faa, "w") as w:
            for g in ids:
                c, st, segs = genes[g]
                if c in seqs:
                    w.write(f">{g}\n{protein(seqs[c], segs, st)}\n")
                else:
                    print(f"WARNING: {g}: contig {c} not in {a.fasta}", file=sys.stderr)
    prots = {}
    for l in open(faa):
        if l.startswith(">"):
            name = l[1:].strip()
        else:
            prots[name] = l.strip()
    ids = [g for g in ids if g in prots]
    print(f"{len(ids)} new genes with a protein ({faa})")

    failed = False
    if not a.summary_only:
        query = f"all [filter] NOT txid{a.exclude_taxid}[Organism:exp]"
        batches = [ids[i:i + a.batch] for i in range(0, len(ids), a.batch)]
        for k, bt in enumerate(batches, 1):
            base = os.path.join(a.out_dir, "batches", f"batch_{k:03d}")
            if os.path.exists(base + ".ids"):
                continue
            q = base + ".faa"
            with open(q, "w") as w:
                for g in bt:
                    w.write(f">{g}\n{prots[g]}\n")
            t0 = time.time()
            cmd = ["blastp", "-remote", "-db", a.db, "-query", q, "-entrez_query", query, "-evalue", a.evalue,
                   "-max_target_seqs", "5", "-outfmt", "6 " + BLAST_COLS, "-out", base + ".part"]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                print(f"batch {k}/{len(batches)} failed ({r.stderr.strip()[:300]}); stopping - submit again to resume",
                      file=sys.stderr)
                failed = True
                break
            os.replace(base + ".part", base + ".tsv")
            open(base + ".ids", "w").write("\n".join(bt) + "\n")   # written last: the batch counts as done
            print(f"batch {k}/{len(batches)} done ({len(bt)} proteins, {time.time() - t0:.0f} s)", flush=True)
            time.sleep(a.pause)
    summary(a.out_dir, ids)
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
