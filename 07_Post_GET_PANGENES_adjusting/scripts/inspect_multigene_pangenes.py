#!/usr/bin/env python3
"""
Pangenes with more than one gene in the same genome: why the genes are together.

get_pangenes groups genes by where they fall on the genome alignments, so two genes of one genome land in the same
pangene when both fall on the place of one gene of the other genomes. For each (pangene, genome) with 2 or more genes,
each gene's protein is aligned (blastp, best HSP) with a reference protein of the pangene - the gene of a genome that
has exactly one gene in the pangene: the reference genome's (--reference) when it has one, otherwise the first genome
of --genomes whose gene has a valid ORF - and the genes are classified:

  tandem_copy   one of the genes is a tandem-array copy added by step 05 (tandem_array_copy=True): expected
  copy          every gene covers >= --min-cov of the reference protein: duplicated copies
  split         the genes cover different parts of the reference protein (they overlap by < 20% of the shorter part)
                and together >= --min-cov of it: one gene annotated in pieces in this genome
  fragment      one gene covers >= --min-cov, another covers less: a piece left beside a full model
  unrelated     a gene has no blastp hit on the reference protein: a different gene lumped into the pangene
  other         none of the above (e.g. pieces that overlap, or all genes partial)
  no_reference  no genome has exactly one gene in the pangene

Outputs (in --out-dir):
  multigene_pangenes.tsv  one row per (pangene, genome): the case, the genes with their basket (step 05 panel gene),
                          orf_status, model_origin, the part of the reference protein each one covers, identity and
                          coverage, and whether the genes are neighbours on the same strand
  summary.tsv             (pangene, genome) rows per case

Needs blastp on the PATH, and step 05's scripts (check_orf.py, decide_models.py) in --step05-scripts.

Usage:
    python3 inspect_multigene_pangenes.py --genes-matrix <pangenes>/pangene_matrix_genes.tr.tab \\
        --gff-dir basket_annotation/gff_prefer_panel_model_final --genome-dir genomes_nonredundant \\
        --genomes genome_order.txt --reference ME49 \\
        --step05-scripts 05_Basket_annotation_integration/scripts --out-dir multigene_inspection
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

EMPTY = {"", "-", "NA"}
CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
                                               "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def protein(cds):
    s = cds.upper()
    return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3)).replace("*", "")


def read_matrix(path):
    """pangene_matrix_genes.tr.tab -> (genomes, {pangene: {genome: [gene IDs]}})."""
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    genomes = [c.replace("_for_get_pangenes", "") for c in rows[0][1:]]
    out = {}
    for r in rows[1:]:
        if not r or r[0].startswith("chr:"):
            continue
        cells = {}
        for g, c in zip(genomes, r[1:]):
            if c not in EMPTY:
                cells[g] = [x for x in re.split(r"[,;]", c) if x]
        out[r[0]] = cells
    return genomes, out


def gene_info(path):
    """gene ID -> dict(contig, strand, start, basket, orf_status, model_origin, tandem)."""
    info = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "gene":
            a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
            info[a["ID"]] = {"contig": f[0], "strand": f[6], "start": int(f[3]), "basket": a.get("basket", ""),
                             "orf_status": a.get("orf_status", ""), "model_origin": a.get("model_origin", ""),
                             "tandem": a.get("tandem_array_copy") == "True"}
    return info


def blastp_best(query, subject):
    """best HSP of query against subject: (identity, subject start, subject end) or None."""
    with tempfile.TemporaryDirectory() as t:
        open(f"{t}/q.faa", "w").write(f">q\n{query}\n")
        open(f"{t}/s.faa", "w").write(f">s\n{subject}\n")
        out = subprocess.run(["blastp", "-query", f"{t}/q.faa", "-subject", f"{t}/s.faa", "-evalue", "1e-5",
                              "-outfmt", "6 pident sstart send bitscore"], capture_output=True, text=True).stdout
    best = None
    for line in out.splitlines():
        p, s0, s1, bs = line.split("\t")
        if best is None or float(bs) > best[3]:
            best = (float(p), int(s0), int(s1), float(bs))
    return best[:3] if best else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--genomes", required=True, help="genome list (order used to pick the reference protein)")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--step05-scripts", required=True, help="folder with step 05's check_orf.py and decide_models.py")
    ap.add_argument("--min-cov", type=float, default=0.7)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    sys.path.insert(0, a.step05_scripts)
    import check_orf                                  # noqa: E402
    from decide_models import read_models, spliced    # noqa: E402

    order = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    _, pg = read_matrix(a.genes_matrix)
    multi = {p: cells for p, cells in pg.items() if any(len(v) > 1 for v in cells.values())}
    print(f"{len(pg)} pangenes, {len(multi)} with more than one gene in a genome", file=sys.stderr)

    # genes whose protein is needed: every gene of the multi cells, and the candidate reference genes
    want = defaultdict(set)
    for p, cells in multi.items():
        for g, ids in cells.items():
            want[g].update(ids)
    info, prot, pos = {}, {}, {}
    for g in order:
        path = os.path.join(a.gff_dir, g + a.gff_suffix)
        gi = gene_info(path)
        info.update({k: dict(v, genome=g) for k, v in gi.items()})
        byc = defaultdict(list)
        for k, v in gi.items():
            byc[(v["contig"], v["strand"])].append((v["start"], k))
        for k2 in byc:
            for i, (_, k) in enumerate(sorted(byc[k2])):
                pos[k] = (k2, i)
        if want[g]:
            genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
            for k, m in read_models(path, want[g]).items():
                prot[k] = protein(spliced(genome, m))
            del genome

    rows, cases = [], Counter()
    for p, cells in sorted(multi.items()):
        ref_gene = None
        singles = {g: ids[0] for g, ids in cells.items() if len(ids) == 1}
        if a.reference in singles:
            ref_gene = singles[a.reference]
        else:
            for g in order:
                if g in singles and info.get(singles[g], {}).get("orf_status") == "valid":
                    ref_gene = singles[g]
                    break
        if ref_gene and ref_gene not in prot:
            g = info[ref_gene]["genome"]
            genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
            for k, m in read_models(os.path.join(a.gff_dir, g + a.gff_suffix), {ref_gene}).items():
                prot[k] = protein(spliced(genome, m))
            del genome
        rp = prot.get(ref_gene, "")
        for g, ids in sorted(cells.items()):
            if len(ids) < 2:
                continue
            hits = {k: (blastp_best(prot[k], rp) if rp and prot.get(k) else None) for k in ids}
            covs = {k: ((h[2] - h[1] + 1) / len(rp) if h else 0.0) for k, h in hits.items()}
            if any(info[k]["tandem"] for k in ids):
                case = "tandem_copy"
            elif not rp:
                case = "no_reference"
            elif any(h is None for h in hits.values()):
                case = "unrelated"
            elif all(c >= a.min_cov for c in covs.values()):
                case = "copy"
            else:
                iv = sorted((h[1], h[2]) for h in hits.values())
                disjoint = all(min(x[1], y[1]) - max(x[0], y[0]) + 1 < 0.2 * min(x[1] - x[0] + 1, y[1] - y[0] + 1)
                               for i, x in enumerate(iv) for y in iv[i + 1:])
                covered = set()
                for s0, s1 in iv:
                    covered.update(range(s0, s1 + 1))
                if disjoint and len(covered) / len(rp) >= a.min_cov:
                    case = "split"
                elif any(c >= a.min_cov for c in covs.values()):
                    case = "fragment"
                else:
                    case = "other"
            cases[case] += 1
            ps = sorted(pos[k] for k in ids)
            neighbours = all(ps[i][0] == ps[i + 1][0] and ps[i + 1][1] - ps[i][1] == 1 for i in range(len(ps) - 1))
            rows.append([p, g, case, ref_gene or "", len(rp), "yes" if neighbours else "no",
                         ";".join(ids),
                         ";".join(info[k]["basket"] for k in ids),
                         ";".join(info[k]["orf_status"] for k in ids),
                         ";".join(info[k]["model_origin"] for k in ids),
                         ";".join(f"{hits[k][1]}-{hits[k][2]}" if hits[k] else "no_hit" for k in ids),
                         ";".join(f"{hits[k][0]:.1f}" if hits[k] else "" for k in ids),
                         ";".join(f"{covs[k]:.2f}" for k in ids)])

    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "multigene_pangenes.tsv"), "w") as fh:
        fh.write("pangene\tgenome\tcase\treference_gene\treference_aa\tneighbours_same_strand\tgenes\tbaskets\t"
                 "orf_status\tmodel_origin\tpart_of_reference_aa\tidentity\tcoverage_of_reference\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        fh.write("case\tpangene_genome_rows\tpangenes\n")
        for c, n in cases.most_common():
            fh.write(f"{c}\t{n}\t{len({r[0] for r in rows if r[2] == c})}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read())


if __name__ == "__main__":
    main()
