#!/usr/bin/env python3
"""
The source-annotation layer of the explorer (e.g. ToxoDB, FungiDB): for every genome, how each gene of the explorer
relates to the genome's own annotation, and where the own genes that left the pangenome were.

For every genome, writes <out-dir>/source_<genome>.json:
  "status":  {final gene: [group, [own gene IDs], model_from, 1 if model_from is a gene of --reference else 0]}
             group: identical (1 own gene, same CDS segment lengths in order) / changed (1 own gene, other CDS) /
                    joined (2+ own genes, now one model) / added (no own gene)
  "orig":    {final gene: [[own ID, start, end, strand, [cds_start, cds_end, ...]], ...]}
  "dropped": {contig: [[own ID, start, end, strand, [cds...], reason], ...]}     own genes with no final gene
  "gap":     [final genes whose ORF is broken by an assembly gap (orf_status=invalid_assembly_gap)]
  "orf":     {own ID: 0 broken ORF / -1 broken by an assembly gap (N)}           with --genome-dir; complete ones left out
Pairing own -> final gene comes from step 05's gene fate table ("kept" -> "<final gene> (<own gene>)"); the reason an
own gene has no final gene is its fate and detail there; with --removed (remove_panel_genes.py's removed_genes.tsv),
an own gene whose final gene was removed in step 07 gets that reason instead (e.g. curation, TE).

Usage:
    python3 build_source_layer.py --gene-fate basket_annotation/gene_fate_prefer_panel_model.tsv \\
        --own-gff-dir basket_annotation/input/native --gff-dir basket_annotation/gff_prefer_panel_model_final \\
        [--genome-dir genomes] [--removed pangenome_final/removed_genes.tsv] -o pangenome_explorer/final/data
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import os
import re

ID_RE, PAR_RE = re.compile(r"(?:^|;)ID=([^;\n]+)"), re.compile(r"(?:^|;)Parent=([^;\n]+)")
att = lambda s, k: (re.search(r"(?:^|;)" + k + r"=([^;\n]*)", s) or [None, ""])[1]
openf = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)
STOP = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def models(path):
    """top-level genes with a CDS: {ID: [ID, start, end, strand, contig, [cds flat]]}"""
    genes, top, cds = {}, {}, collections.defaultdict(set)
    for l in openf(path):
        if l.startswith("##FASTA"):
            break
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        i, p = ID_RE.search(f[8]), PAR_RE.search(f[8])
        if not p:
            if i:
                genes[i.group(1)] = [i.group(1), int(f[3]), int(f[4]), f[6], f[0]]
                top[i.group(1)] = i.group(1)
            continue
        g = top.get(p.group(1).split(",")[0])
        if g is None:
            continue
        if i:
            top[i.group(1)] = g
        if f[2] == "CDS":
            cds[g].add((int(f[3]), int(f[4])))
    for g, x in genes.items():
        x.append([v for s, e in sorted(cds[g]) for v in (s, e)])
    return {g: x for g, x in genes.items() if x[5]}


def seglens(flat):
    return [flat[i + 1] - flat[i] for i in range(0, len(flat), 2)]


def fasta(path):
    seqs, name, parts = {}, None, []
    for l in openf(path):
        if l.startswith(">"):
            if name:
                seqs[name] = "".join(parts)
            name, parts = l[1:].split()[0], []
        else:
            parts.append(l.strip())
    if name:
        seqs[name] = "".join(parts)
    return seqs


def orf(seq, flat, strand):
    """1 complete, 0 broken, -1 broken by an assembly gap"""
    s = "".join(seq[flat[i] - 1:flat[i + 1]] for i in range(0, len(flat), 2)).upper()
    if strand == "-":
        s = s.translate(COMP)[::-1]
    if "N" in s:
        return -1
    if len(s) % 3 or not s.startswith("ATG") or s[-3:] not in STOP:
        return 0
    return 0 if any(s[i:i + 3] in STOP for i in range(0, len(s) - 3, 3)) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-fate", required=True, help="step 05's gene_fate_<mode>.tsv")
    ap.add_argument("--own-gff-dir", required=True, help="each genome's own annotation as given to step 05 (<genome>.gff3/.gff)")
    ap.add_argument("--gff-dir", required=True, help="the GFFs of the explorer (<genome>_basket_annotation.gff)")
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--genome-dir", default="", help="<genome>.fna: the ORF of each own gene (optional)")
    ap.add_argument("--reference", default="", help="the reference genome (its own genes mark 'added from the reference')")
    ap.add_argument("--removed", default="", help="remove_panel_genes.py's removed_genes.tsv (final set; optional)")
    ap.add_argument("-o", "--out-dir", required=True)
    a = ap.parse_args()

    fate = collections.defaultdict(dict)
    for r in csv.DictReader(open(a.gene_fate), delimiter="\t"):
        fate[r["genome"]][r["gene"]] = (r["fate"], r["detail"])
    removed = {}                                   # (genome, final gene) -> reason it was removed in step 07
    if a.removed:
        for r in csv.DictReader(open(a.removed), delimiter="\t"):
            removed[(r["genome"], r["gene"])] = r["reason"]
    os.makedirs(a.out_dir, exist_ok=True)
    ref_ids = set(fate.get(a.reference, {}))
    for fgff in sorted(glob.glob(os.path.join(a.gff_dir, "*" + a.gff_suffix))):
        G = os.path.basename(fgff)[:-len(a.gff_suffix)]
        own_files = sorted(glob.glob(os.path.join(a.own_gff_dir, G + ".gff*")))
        if not own_files or G not in fate:
            print(f"{G}: no own annotation or gene fate - skipped")
            continue
        own = models(own_files[0])
        fin, fattr = {}, {}
        for l in openf(fgff):
            f = l.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] == "gene":
                fattr[att(f[8], "ID")] = (att(f[8], "model_from"), att(f[8], "orf_status"))
        fin = models(fgff)
        by_final = collections.defaultdict(list)
        for g, (fa, det) in fate[G].items():
            if fa == "kept" and det:
                by_final[det.split()[0]].append(g)
        status, orig, used = {}, {}, set()
        for fg, (mfrom, _) in fattr.items():
            gs = [g for g in by_final.get(fg, []) if g in own]
            used.update(gs)
            if not gs:
                grp = "added"
            elif len(gs) > 1:
                grp = "joined"
            else:
                grp = "identical" if fg in fin and seglens(fin[fg][5]) == seglens(own[gs[0]][5]) else "changed"
            status[fg] = [grp, gs, mfrom, 1 if mfrom in ref_ids else 0]
            if gs:
                orig[fg] = [[own[g][0], own[g][1], own[g][2], own[g][3], own[g][5]] for g in gs]
        dropped = collections.defaultdict(list)
        for g, x in own.items():
            if g in used:
                continue
            fa, det = fate[G].get(g, ("not in the pangenome", ""))
            why = f"{fa}: {det}".strip(": ")
            if fa == "kept" and det and (G, det.split()[0]) in removed:
                why = "removed in step 07 - " + removed[(G, det.split()[0])]
            dropped[x[4]].append([x[0], x[1], x[2], x[3], x[5], why])
        for L in dropped.values():
            L.sort(key=lambda r: r[1])
        gap = [g for g, (_, o) in fattr.items() if o == "invalid_assembly_gap"]
        torf = {}
        if a.genome_dir:
            fa_files = sorted(glob.glob(os.path.join(a.genome_dir, G + ".f*a*")))
            if fa_files:
                seqs = fasta(fa_files[0])
                for g, x in own.items():
                    if x[4] in seqs:
                        v = orf(seqs[x[4]], x[5], x[3])
                        if v != 1:
                            torf[g] = v
        out = os.path.join(a.out_dir, f"source_{G}.json")
        json.dump({"status": status, "orig": orig, "dropped": dropped, "gap": gap, "orf": torf}, open(out, "w"),
                  separators=(",", ":"))
        n = collections.Counter(v[0] for v in status.values())
        print(f"{G}: {dict(n)}; own genes not in the pangenome {sum(len(v) for v in dropped.values())}; "
              f"own models with a broken ORF {len(torf) if a.genome_dir else 'not checked'}")


if __name__ == "__main__":
    main()
