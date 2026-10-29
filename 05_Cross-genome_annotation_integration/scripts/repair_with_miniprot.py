#!/usr/bin/env python3
"""
Try to repair, with miniprot, the genes written with an invalid ORF by build_basket_annotation.py (orf_invalid=True),
except selenoproteins and pseudogenes (pseudogene=True, from prepare_pseudogenes.py: orf_status=pseudogene).

For each such gene, the protein of its panel gene (attribute basket=: the panel gene's own model, in its own genome) is
aligned by miniprot to the genome the gene is in. The miniprot model that overlaps the gene's place (same contig and
strand, >= 50% of the gene) is taken; its CDS are checked (ATG start, stop, length x3, no internal stop):
  repaired          valid ORF and >= --min-protein-cov of the protein aligned: the gene's model is replaced by
                    miniprot's (same gene ID; model_origin=repaired_miniprot; orf_invalid removed)
  frameshift        miniprot needed a frameshift (fs > 0): flag frameshift=<count> (an indel, often an assembly error)
  stop_codon        miniprot crossed an in-frame stop (st > 0): flag stop_codon=<count> (pseudogene, or a sequencing error)
  start/stop_codon_lost  the whole protein aligns, but the genome has no ATG at the start or no stop at the end
  partial           the protein is not covered to the end; flag partial_alignment=True
  repair_overlaps_another_gene  a valid model was found but it would lie over another gene: not used
  not_tried_assembly_gap / not_tried_near_contig_end   not sent to miniprot: the gene has an N in its CDS, or lies
                    within 1 kb of a contig end (the assembly explains the broken ORF): orf_status=invalid_assembly_gap /
                    invalid_near_contig_end
Every gene line gets orf_status=valid / repaired_miniprot / selenoprotein / pseudogene / invalid_assembly_gap /
invalid_near_contig_end / invalid (one value per gene, to filter or colour by).
  repair_structure_differs      a valid model was found, but its exons are not the panel model's (not the same number,
                    or an exon differs by more than --max-exon-diff bp): not used, the gene keeps its transferred model
  repair_noncanonical_intron    a valid model was found but one of its introns has no canonical splice sites (GT/GC..AG):
                    not used (miniprot spliced out an insertion, from the assembly or a real mutation)
  no_hit            miniprot found nothing at the gene's place
Writes the GFFs again in --out-dir (same names), with miniprot_results.tsv (one row per gene tried).

Usage: python3 repair_with_miniprot.py --gff-dir gff_similar --native-dir gff_for_liftoff --genome-dir genomes \\
           --genomes genome_order.txt --out-dir gff_similar_miniprot [--miniprot miniprot] [--threads 8]
Needs miniprot on the PATH (or --miniprot).
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_orf   # noqa: E402
from decide_models import read_models, spliced   # noqa: E402

BASES = "TCAG"
# Liftoff's ORF marks describe the model miniprot replaces: dropped from a repaired gene (every other attribute is kept)
STALE_ORF_TAGS = ("valid_ORF", "valid_ORFs", "missing_start_codon", "missing_stop_codon", "inframe_stop_codon")
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES)
         for k, c in enumerate(BASES)}


def with_orf_status(col9, why=""):
    """add orf_status= (one value per gene, easy to filter or colour by); why = assembly_gap / near_contig_end for a gene
    not sent to miniprot because the assembly explains its broken ORF."""
    if "orf_invalid=True" not in col9:
        status = "valid"
    elif "selenoprotein=True" in col9:
        status = "selenoprotein"
    elif "pseudogene=True" in col9:
        status = "pseudogene"
    elif why:
        status = "invalid_" + why
    else:
        status = "invalid"
    return col9 + f";orf_status={status}"


def exon_lengths(segs, strand):
    """CDS exon lengths in coding order."""
    L = [y - x + 1 for x, y in sorted(segs)]
    return L[::-1] if strand == "-" else L


def same_structure(rep, panel, max_diff):
    """the repaired model has the panel model's exons: same number, each within max_diff bp."""
    return bool(panel) and len(rep) == len(panel) and all(abs(x - y) <= max_diff for x, y in zip(rep, panel))


def canonical_intron(seq, s, e, strand):
    """the intron s..e (1-based, inclusive) starts with GT or GC and ends with AG, on the gene's strand."""
    x = seq[s - 1:e].upper()
    if strand == "-":
        x = x.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    return len(x) >= 4 and x[:2] in ("GT", "GC") and x[-2:] == "AG"


def translate(seq):
    p = "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))
    return p[:-1] if p.endswith("*") else p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--suffix", default="_basket_annotation.gff")
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--miniprot", default="miniprot")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-protein-cov", type=float, default=0.9)
    ap.add_argument("--max-exon-diff", type=int, default=10,
                    help="a repair keeps the panel model's structure: same exons, each within this many bp")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]

    # genes to try, per genome
    todo = defaultdict(dict)   # genome -> gene ID -> (contig, strand, start, end, basket)
    own_cds = defaultdict(list)   # (genome, gene ID) -> CDS segments as written
    others = defaultdict(list)    # (genome, contig, strand) -> (CDS start, CDS end, gene ID) of every gene
    for g in genomes:
        cds_of = defaultdict(list)
        for line in open(os.path.join(a.gff_dir, g + a.suffix)):
            f = line.rstrip("\n").split("\t")
            if (len(f) > 8 and f[2] == "gene" and "orf_invalid=True" in f[8] and "selenoprotein=True" not in f[8]
                    and "pseudogene=True" not in f[8]):
                gid = re.search(r"ID=([^;]+)", f[8]).group(1)
                todo[g][gid] = (f[0], f[6], int(f[3]), int(f[4]), re.search(r"basket=([^;]+)", f[8]).group(1))
            if len(f) > 8 and f[2] == "CDS":
                gid = re.search(r"Parent=([^;]+)", f[8]).group(1).split(".t")[0]
                cds_of[gid].append((int(f[3]), int(f[4]), f[0], f[6]))
        for gid, c in cds_of.items():
            own_cds[(g, gid)] = [(x, y) for x, y, _, _ in c]
            others[(g, c[0][2], c[0][3])].append((min(x for x, *_ in c), max(y for _, y, *_ in c), gid))
    panel = {b for g in todo for *_, b in todo[g].values()}
    # proteins of the panel genes
    prot, owner, panel_exons = {}, {}, {}
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))
        models = read_models(path, panel)
        if not models:
            continue
        genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
        for k, m in models.items():
            prot[k] = translate(spliced(genome, m)).replace("*", "X")
            owner[k] = g
            panel_exons[k] = exon_lengths(m[2], m[1])
    results, new_models = [], {}
    for g in genomes:
        if not todo[g]:
            continue
        genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
        # genes broken by the assembly are not sent to miniprot: an N in the CDS, or the gene near a contig end
        skip = {}
        for gid, (c, st, s, e, b) in todo[g].items():
            seq = genome.get(c, "")
            if any("N" in seq[x - 1:y].upper() for x, y in own_cds.get((g, gid), [])):
                skip[gid] = "assembly_gap"
            elif s < 1000 or e > len(seq) - 1000:
                skip[gid] = "near_contig_end"
        for gid, why in skip.items():
            results.append([g, gid, todo[g][gid][4], "not_tried_" + why, "", "", "", why])
            del todo[g][gid]
        with tempfile.TemporaryDirectory() as tmp:
            faa = os.path.join(tmp, "q.faa")
            with open(faa, "w") as fh:
                for gid, (_, _, _, _, b) in todo[g].items():
                    if prot.get(b):
                        fh.write(f">{gid}\n{prot[b]}\n")
            out = subprocess.run([a.miniprot, "-t", str(a.threads), "--gff", "-N", "5",
                                  os.path.join(a.genome_dir, g + a.genome_suffix), faa],
                                 check=True, capture_output=True, text=True).stdout
        hits = defaultdict(list)   # gene ID -> [(paf tags, mRNA fields, CDS segments)]
        paf, cur = None, None
        for line in out.splitlines():
            if line.startswith("##PAF"):
                paf = line.split("\t")
                continue
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                tags = dict((t.split(":")[0], t.split(":")[-1]) for t in paf[13:] if t.count(":") >= 2)
                tgt = re.search(r"Target=(\S+) (\d+) (\d+)", f[8])
                cur = {"q": tgt.group(1), "qs": int(tgt.group(2)), "qe": int(tgt.group(3)), "qlen": int(paf[2]),
                       "fs": int(tags.get("fs", 0)), "st": int(tags.get("st", 0)), "c": f[0], "strand": f[6],
                       "cds": []}
                hits[cur["q"]].append(cur)
            elif f[2] in ("CDS", "stop_codon") and cur:
                cur["cds"].append((int(f[3]), int(f[4])))
        for gid, (c, st, s, e, b) in todo[g].items():
            best = None
            for h in hits.get(gid, []):
                if h["c"] != c or h["strand"] != st or not h["cds"]:
                    continue
                hs, he = min(x for x, _ in h["cds"]), max(y for _, y in h["cds"])
                if min(e, he) - max(s, hs) + 1 >= 0.5 * (e - s + 1):
                    if best is None or (h["qe"] - h["qs"]) > (best["qe"] - best["qs"]):
                        best = h
            row = [g, gid, b]
            if not best:
                results.append(row + ["no_hit", "", "", "", ""])
                continue
            pcov = (best["qe"] - best["qs"] + 1) / best["qlen"]
            segs = []
            for x, y in sorted(best["cds"]):   # join the stop codon to its CDS when they touch
                if segs and x <= segs[-1][1] + 1:
                    segs[-1] = (segs[-1][0], max(segs[-1][1], y))
                else:
                    segs.append((x, y))
            start_ok, stop_ok, frame_ok, internal, ok = check_orf.check(spliced(genome, (c, st, segs)))
            seq = genome.get(c, "")
            gap = any("N" in seq[x - 1:y] for x, y in own_cds.get((g, gid), []))
            ns, ne = segs[0][0], segs[-1][1]
            clash = next((o for x, y, o in others.get((g, c, st), [])
                          if o != gid and min(ne, y) - max(ns, x) + 1 >= 0.5 * min(ne - ns + 1, y - x + 1)), None)
            near_end = s < 1000 or e > len(seq) - 1000
            noncanon = sum(1 for (x1, y1), (x2, y2) in zip(segs, segs[1:]) if not canonical_intron(seq, y1 + 1, x2 - 1, st))
            if ok and pcov >= a.min_protein_cov and clash:
                cls = "repair_overlaps_another_gene"
            elif ok and pcov >= a.min_protein_cov and not same_structure(exon_lengths(segs, st), panel_exons.get(b),
                                                                         a.max_exon_diff):
                cls = "repair_structure_differs"   # a valid ORF, but not the panel model's exon structure: not used
            elif ok and pcov >= a.min_protein_cov and noncanon:
                cls = "repair_noncanonical_intron"   # miniprot spliced out sequence with no splice signal: not used
            elif ok and pcov >= a.min_protein_cov:
                cls = "repaired"
                new_models[(g, gid)] = (c, st, segs)
            elif best["fs"] > 0:
                cls = "frameshift"
            elif best["st"] > 0:
                cls = "stop_codon"
            elif pcov < a.min_protein_cov:
                cls = "partial"
            elif not start_ok and not stop_ok:
                cls = "start_and_stop_codon_lost"
            elif not start_ok:
                cls = "start_codon_lost"
            elif not stop_ok:
                cls = "stop_codon_lost"
            else:
                cls = "invalid_orf_other"
            results.append(row + [cls, f"{pcov:.3f}", best["fs"], best["st"],
                                  ";".join(x for x, y in (("assembly_gap", gap), ("near_contig_end", near_end)) if y)])
        print(f"  {g}: {len(todo[g])} genes tried, {len(skip)} not sent (N in the CDS or near a contig end)", file=sys.stderr)

    # write the GFFs again
    flags = {(r[0], r[1]): r for r in results}
    for g in genomes:
        src = os.path.join(a.gff_dir, g + a.suffix)
        lines = open(src).read().splitlines()
        # attributes of the old mRNA and CDS lines of each gene (all but ID/Parent), carried to a repaired model
        old_attrs = {}
        for line in lines:
            f = line.split("\t")
            if len(f) > 8 and f[2] in ("mRNA", "CDS"):
                par = re.search(r"Parent=([^;]+)", f[8])
                if par:
                    gid0 = par.group(1).split(".t")[0]
                    kept = ";".join(kv for kv in f[8].split(";") if kv and kv.split("=", 1)[0] not in ("ID", "Parent") + STALE_ORF_TAGS)
                    old_attrs.setdefault((gid0, f[2]), kept)
        out_lines, skip = [], None
        for line in lines:
            f = line.split("\t")
            if len(f) < 9:
                out_lines.append(line)
                continue
            if f[2] == "gene":
                gid = re.search(r"ID=([^;]+)", f[8]).group(1)
                skip = None
                r = flags.get((g, gid))
                if r and (g, gid) in new_models:
                    c, st, segs = new_models[(g, gid)]
                    f[3], f[4] = str(segs[0][0]), str(segs[-1][1])
                    attrs = ";".join(kv for kv in f[8].split(";") if kv.split("=", 1)[0] not in STALE_ORF_TAGS + ("orf_invalid",))
                    attrs = re.sub(r"model_origin=[^;]+", "model_origin=repaired_miniprot", attrs)
                    f[8] = attrs + ";orf_status=repaired_miniprot"
                    out_lines.append("\t".join(f))
                    tid = f"{gid}.t1"
                    def with_old(base, kind):
                        return ";".join(x for x in (base, old_attrs.get((gid, kind), "")) if x)
                    out_lines.append("\t".join(f[:2] + ["mRNA"] + f[3:8] + [with_old(f"ID={tid};Parent={gid}", "mRNA")]))
                    for n, (x, y) in enumerate(segs, 1):
                        out_lines.append("\t".join([c, f[1], "exon", str(x), str(y), ".", st, ".",
                                                    f"ID={tid}.exon{n};Parent={tid}"]))
                    for n, (x, y) in enumerate(segs, 1):
                        out_lines.append("\t".join([c, f[1], "CDS", str(x), str(y), ".", st, ".",
                                                    with_old(f"ID={tid}.cds{n};Parent={tid}", "CDS")]))
                    skip = gid
                    continue
                if r:
                    extra = {"frameshift": f"frameshift={r[5]}", "stop_codon": f"stop_codon={r[6]}",
                             "repair_overlaps_another_gene": "miniprot_repair_overlaps_another_gene=True",
                             "repair_noncanonical_intron": "miniprot_repair_noncanonical_intron=True",
                             "repair_structure_differs": "miniprot_repair_structure_differs=True",
                             "partial": "partial_alignment=True", "no_hit": "no_miniprot_hit=True",
                             "start_codon_lost": "start_codon_lost=True", "stop_codon_lost": "stop_codon_lost=True",
                             "start_and_stop_codon_lost": "start_codon_lost=True;stop_codon_lost=True"}.get(r[3], "")
                    f[8] = ";".join(x for x in (f[8], extra) if x)
                f[8] = with_orf_status(f[8], r[7] if r else "")
                out_lines.append("\t".join(f))
            else:
                if skip and re.search(r"Parent=" + re.escape(skip) + r"(\.|;|$)", f[8]):
                    continue
                out_lines.append(line)
        with open(os.path.join(a.out_dir, g + a.suffix), "w") as fh:
            fh.write("\n".join(out_lines) + "\n")
    with open(os.path.join(a.out_dir, "miniprot_results.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "gene", "panel_gene", "result", "protein_cov", "frameshifts", "stop_codons", "context"])
        w.writerows(results)
    from collections import Counter
    c = Counter(r[3] for r in results)
    print(f"genes with an invalid ORF: {len(results)} (sent to miniprot: {sum(1 for r in results if not r[3].startswith('not_tried'))})")
    for k, v in c.most_common():
        print(f"  {v:6}  {k}")


if __name__ == "__main__":
    main()
