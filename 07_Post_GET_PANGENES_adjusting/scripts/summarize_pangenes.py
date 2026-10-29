#!/usr/bin/env python3
"""
Two lists to help decide which genes to remove after curation (nothing is removed here):

  pangene_quality.tsv          one row per pangene: in how many strains its genes have an invalid ORF, where the ORF
                               breaks and whether the strains share the same defect, the PSAURON score of its panel
                               gene(s) and its InterPro domains; sorted by the share of strains with an invalid ORF.
                               A strain counts as valid when at least one of its genes in the pangene has a valid ORF
                               (valid, repaired_miniprot, selenoprotein), invalid when its genes are all invalid or
                               pseudogenes, and is not judged when its ORFs are broken by the assembly (N in the CDS,
                               contig end). Pangenes with an invalid ORF in more than --min-share of the judged strains
                               are marked low_confidence=yes, or low_confidence=pseudogene when they hold a reference
                               pseudogene (a known pseudogene).
                               The defect of each invalid gene (with --genome-dir) is read on its CDS: start_lost (no
                               ATG), stop_lost (no stop at the end), internal_stop@aaN (first stop codon inside, at codon
                               N), frameshift(len%3=R) (CDS length not a multiple of 3). The genes of a panel gene come
                               from the same model, so the same defect in several strains is likely one shared mutation
                               (or an error of the model); different defects point to independent losses.
  genes_not_in_pangenome.tsv   genes of the GFFs that are in no pangene (get_pangenes leaves out CDS < 100 bp), with
                               their panel gene, CDS length, ORF status and PSAURON.
  summary.tsv                  counts

Inputs: the pangene matrix (pangene_matrix_genes.tr.tab), the GFFs given to get_pangenes (step 05: basket=, orf_status,
model_origin, Dbxref InterPro), step 05's winners.tsv (PSAURON of every panel gene), and the genomes (--genome-dir,
<genome>.fna; optional, for the defects) and InterPro entry.list (--entry-list; optional, for domain names).

Usage:
    python3 summarize_pangenes.py --genes-matrix <pangenes>/pangene_matrix_genes.tr.tab \\
        --gff-dir basket_annotation/gff_prefer_panel_model_final --winners basket_annotation/winners.tsv \\
        --genome-dir genomes --entry-list interpro_entry.list --out-dir pangenome_fixed/pangene_lists
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

from fix_pangenome_safe import read_matrix

VALID = {"valid", "repaired_miniprot", "selenoprotein"}
INVALID = {"invalid", "pseudogene"}
ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def attr(a, k):
    m = re.search(rf"(?:^|;){k}=([^;\n]*)", a)
    return m.group(1) if m else ""


def read_genes(gff_dir, suffix):
    """gene ID -> dict(genome, location, basket, orf_status, model_origin, description, cds_bp, strand, cds (segments
    of its first transcript), ipr (InterPro IDs))."""
    genes = {}
    for path in sorted(glob.glob(os.path.join(gff_dir, "*" + suffix))):
        genome = os.path.basename(path)[:-len(suffix)]
        mrna, first = {}, {}
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "gene":
                gid = ID_RE.search(f[8]).group(1)
                genes[gid] = {"genome": genome, "location": f"{f[0]}:{f[3]}-{f[4]}{f[6]}",
                              "basket": attr(f[8], "basket"), "orf_status": attr(f[8], "orf_status"),
                              "model_origin": attr(f[8], "model_origin"),
                              "description": attr(f[8], "description").replace("%2C", ","), "cds_bp": 0,
                              "strand": f[6], "cds": [], "ipr": set()}
            elif f[2] in ("mRNA", "transcript"):
                i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
                if i and p:
                    mrna[i.group(1)] = p.group(1)
                    first.setdefault(p.group(1), i.group(1))
            elif f[2] == "CDS":
                p = PARENT_RE.search(f[8])
                g = mrna.get(p.group(1), p.group(1)) if p else None
                if g in genes and first.get(g, p.group(1)) == p.group(1):
                    genes[g]["cds_bp"] += int(f[4]) - int(f[3]) + 1
                    genes[g]["cds"].append((f[0], int(f[3]), int(f[4])))
                    d = attr(f[8], "Dbxref")
                    genes[g]["ipr"].update(x.split(":", 1)[1] for x in d.split(",") if x.startswith("InterPro:"))
    return genes


def read_contigs(path, wanted):
    """the sequences of the contigs in wanted, from a FASTA file."""
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name in wanted:
                seqs[name] = "".join(buf).upper()
            name, buf = line[1:].split()[0], []
        elif name in wanted:
            buf.append(line.strip())
    if name in wanted:
        seqs[name] = "".join(buf).upper()
    return seqs


def defect(cds):
    """where the ORF of a CDS breaks: start_lost, internal_stop@aaN, frameshift(len%3=R), stop_lost; or none."""
    out = []
    if not cds.startswith("ATG"):
        out.append("start_lost")
    codons = [cds[i:i + 3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    inside = [i for i, c in enumerate(codons[:-1]) if c in STOPS]
    if inside:
        out.append(f"internal_stop@aa{inside[0] + 1}")
    if len(cds) % 3:
        out.append(f"frameshift(len%3={len(cds) % 3})")
    elif not codons or codons[-1] not in STOPS:
        out.append("stop_lost")
    return ";".join(out) or "none"


def defects(genes, genome_dir, suffix):
    """gene ID -> defect, for the genes with orf_status=invalid."""
    want = defaultdict(list)
    for g, x in genes.items():
        if x["orf_status"] == "invalid" and x["cds"]:
            want[x["genome"]].append(g)
    out = {}
    for genome, ids in want.items():
        path = os.path.join(genome_dir, genome + suffix)
        if not os.path.exists(path):
            print(f"WARNING: no {path}; no defects for {genome}", file=sys.stderr)
            continue
        seqs = read_contigs(path, {c for g in ids for c, _, _ in genes[g]["cds"]})
        for g in ids:
            x = genes[g]
            segs = sorted(x["cds"], key=lambda t: t[1])
            cds = "".join(seqs.get(c, "")[s - 1:e] for c, s, e in segs)
            if x["strand"] == "-":
                cds = cds.translate(COMP)[::-1]
            out[g] = defect(cds)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--winners", required=True, help="step 05's winners.tsv (PSAURON of every panel gene)")
    ap.add_argument("--genome-dir", default="", help="<genome>.fna, to read where each invalid ORF breaks (optional)")
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--entry-list", default="", help="InterPro entry.list, for the domain names (optional)")
    ap.add_argument("--min-share", type=float, default=0.5,
                    help="low confidence: invalid ORF in more than this share of the judged strains")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    genes = read_genes(a.gff_dir, a.gff_suffix)
    dfc = defects(genes, a.genome_dir, a.genome_suffix) if a.genome_dir else {}
    names = {}
    if a.entry_list and os.path.exists(a.entry_list):
        for line in open(a.entry_list):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0].startswith("IPR"):
                names[f[0]] = f[2]
    psauron = {r["gene"]: r.get("psauron", "") for r in csv.DictReader(open(a.winners), delimiter="\t")}
    header, genomes, order, cells, extra = read_matrix(a.genes_matrix)
    os.makedirs(a.out_dir, exist_ok=True)

    # 1. pangene quality
    rows, in_matrix = [], set()
    for p in order:
        status = {}
        baskets, desc, ipr = Counter(), Counter(), Counter()
        pseudo = False
        strain_defect = {}
        for g, ids in cells[p].items():
            st = [genes[x]["orf_status"] for x in ids if x in genes]
            pseudo = pseudo or "pseudogene" in st
            in_matrix.update(ids)
            for x in ids:
                if x in genes:
                    baskets[genes[x]["basket"]] += 1
                    if genes[x]["description"]:
                        desc[genes[x]["description"]] += 1
                    ipr.update(genes[x]["ipr"])
            status[g] = ("valid" if any(s in VALID for s in st) else
                         "invalid" if st and all(s in INVALID for s in st) else "not_judged")
            if status[g] == "invalid":
                d = [dfc[x] for x in ids if x in dfc]
                if d:
                    strain_defect[g] = d[0]
        c = Counter(status.values())
        judged = c["valid"] + c["invalid"]
        share = c["invalid"] / judged if judged else 0.0
        main = baskets.most_common(1)[0][0] if baskets else ""
        ps = psauron.get(main, "")
        dc = Counter(strain_defect.values())
        top, ntop = dc.most_common(1)[0] if dc else ("", 0)
        n_genes = sum(len(v) for v in cells[p].values())
        rows.append([p, len(cells[p]), n_genes, c["valid"], c["invalid"],
                     c["not_judged"], f"{share:.2f}", ("pseudogene" if pseudo else "yes") if judged and share > a.min_share else "no",
                     top, f"{ntop} of {len(strain_defect)}" if dc else "", len(dc) if dc else "",
                     main, ps, "yes" if ps and float(ps) < 0.5 else "no",
                     "yes" if ipr else "no",
                     ";".join(f"{i} {names.get(i, '')}".strip() + f" ({n}/{n_genes})" for i, n in ipr.most_common()),
                     ";".join(f"{b}:{psauron.get(b, '')}:{n}" for b, n in baskets.most_common()),
                     desc.most_common(1)[0][0] if desc else "",
                     ",".join(sorted(g.replace("_for_get_pangenes", "") for g, s in status.items() if s == "invalid"))])
    rows.sort(key=lambda r: (-float(r[6]), r[0]))
    with open(os.path.join(a.out_dir, "pangene_quality.tsv"), "w") as fh:
        fh.write("pangene\tstrains\tgenes\tstrains_valid\tstrains_invalid\tstrains_not_judged\tshare_invalid\t"
                 "low_confidence\tmost_common_defect\tstrains_with_that_defect\tdistinct_defects\tpanel_gene\tpsauron\t"
                 "psauron_below_0.5\thas_domain\tinterpro(genes with it/genes)\tpanel_genes(id:psauron:genes)\t"
                 "description\tinvalid_strains\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    # 2. genes in no pangene
    out = sorted((g for g in genes if g not in in_matrix), key=lambda g: (genes[g]["genome"], g))
    with open(os.path.join(a.out_dir, "genes_not_in_pangenome.tsv"), "w") as fh:
        fh.write("genome\tgene\tpanel_gene\tlocation\tcds_bp\torf_status\tmodel_origin\tpsauron\tdescription\n")
        for g in out:
            x = genes[g]
            fh.write(f"{x['genome']}\t{g}\t{x['basket']}\t{x['location']}\t{x['cds_bp']}\t{x['orf_status']}\t"
                     f"{x['model_origin']}\t{psauron.get(x['basket'], '')}\t{x['description']}\n")

    lf = [r for r in rows if r[7] == "yes"]
    shared = [r for r in lf if r[9] and 2 <= int(r[9].split()[0]) == int(r[9].split()[2])]
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as fh:
        fh.write(f"pangenes\t{len(rows)}\nlow_confidence (invalid ORF in > {a.min_share:.0%} of strains)\t{len(lf)}\n"
                 f"low_confidence_and_psauron_below_0.5\t{sum(1 for r in lf if r[13] == 'yes')}\n"
                 f"low_confidence_without_domain\t{sum(1 for r in lf if r[14] == 'no')}\n"
                 f"low_confidence_same_defect_in_all_invalid_strains (2+)\t{len(shared)}\n"
                 f"low_confidence_pseudogene (reference pseudogenes, not counted above)\t"
                 f"{sum(1 for r in rows if r[7] == 'pseudogene')}\n"
                 f"pangenes_psauron_below_0.5\t{sum(1 for r in rows if r[13] == 'yes')}\n"
                 f"genes_not_in_pangenome\t{len(out)}\n"
                 f"genes_not_in_pangenome_cds_below_100bp\t{sum(1 for g in out if genes[g]['cds_bp'] < 100)}\n")
    print(open(os.path.join(a.out_dir, "summary.tsv")).read())


if __name__ == "__main__":
    main()
