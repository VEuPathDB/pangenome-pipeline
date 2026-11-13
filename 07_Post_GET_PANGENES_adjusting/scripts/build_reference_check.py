#!/usr/bin/env python3
"""
Build reference_check.html: a page to curate the genes the pipeline added to the reference genome
(not_in_reference_annotation=True in --gff), with a choice per gene that is downloaded as curation_actions.tsv.

Evidence, all optional (a section is left out when its input is not given):
  --junctions       intron junctions GFF (VEuPathDB unifiedIntronJunction): drawn as arcs; for a new gene, how many of
                    its introns a junction confirms (both ends within 3 bp)
  --region-rnaseq   region_rnaseq.py's table: RNA-seq over each new gene's CDS, and the level of a typical
                    reference gene (drawn as a dashed line)
  --winners         step 05's winners.tsv: the PSAURON score of each panel gene
  --quality         part 2's pangene_lists/pangene_quality.tsv: in how many strains the pangene's ORF is valid
  --protein-evidence  protein_evidence.py's table: the best blastp hit of each new protein outside the species
  --genome-fasta    the reference genome: each new gene gets a link that opens NCBI blastp (nr) with its protein,
                    the species of --blast-exclude left out; nothing is searched until the BLAST button is pressed
Plain or gzipped (.gz) GFFs. Run by reference_check_job.sh.
"""
import argparse
import csv
import gzip
import json
import re
import urllib.parse
from collections import defaultdict

FLANK = 1500
MAXSIDE = 20000   # bp shown at most on each side of the gene
att = lambda s, k: (re.search(r"(?:^|;)" + k + r"=([^;\n]*)", s) or [None, ""])[1]
dec = lambda v: urllib.parse.unquote(v or "").replace("+", " ")
openf = lambda p: gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_gff(path):
    """gene ID -> dict(c, s, e, st, desc, attrs of a new gene, cds=[(s, e)]); per-contig sorted gene list"""
    genes, t2g, cds = {}, {}, defaultdict(list)
    for l in openf(path):
        if l.startswith("#") or l.count("\t") < 8:
            continue
        f = l.rstrip("\n").split("\t")
        if f[2] == "gene":
            genes[att(f[8], "ID")] = dict(c=f[0], s=int(f[3]), e=int(f[4]), st=f[6], desc=dec(att(f[8], "description")),
                                          new=att(f[8], "not_in_reference_annotation") == "True",
                                          basket=att(f[8], "basket"), model_from=att(f[8], "model_from"),
                                          origin=att(f[8], "model_origin"), orf=att(f[8], "orf_status"),
                                          domain=att(f[8], "has_domain"))
        elif f[2] in ("mRNA", "transcript"):
            t2g[att(f[8], "ID")] = att(f[8], "Parent")
        elif f[2] == "CDS":
            for p in att(f[8], "Parent").split(","):
                cds[p].append((int(f[3]), int(f[4])))
    for t, segs in cds.items():
        g = t2g.get(t, t)
        if g in genes and len(segs) > len(genes[g].get("cds", [])):
            genes[g]["cds"] = sorted(segs)
    by = defaultdict(list)
    for g, x in genes.items():
        x.setdefault("cds", [])
        by[x["c"]].append((x["s"], x["e"], g))
    for c in by:
        by[c].sort()
    return genes, by


def read_matrix(path, genome):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    cols = [i for i, h in enumerate(rows[0]) if i and h.strip()]
    ci = [i for i in cols if rows[0][i].startswith(genome)][0]
    pg, occ = {}, {}
    for r in rows[1:]:
        if not r or r[0].startswith("chr"):
            continue
        occ[r[0]] = sum(1 for i in cols if i < len(r) and r[i] not in ("", "-", "NA"))
        if ci < len(r):
            for x in re.split(r"[,;]", r[ci]):
                if x and x != "NA":
                    pg[x] = r[0]
    return pg, occ


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gff", required=True, help="the reference's final GFF (step 05)")
    ap.add_argument("--matrix", required=True, help="step 07's pangene_matrix_genes.tr.tab")
    ap.add_argument("--winners", default="")
    ap.add_argument("--quality", default="")
    ap.add_argument("--explorer", default="", help="path or URL of the explorer page (explorer_standalone.html); each pangene ID links to it")
    ap.add_argument("--junctions", default="")
    ap.add_argument("--region-rnaseq", default="")
    ap.add_argument("--rnaseq-calls", default="")
    ap.add_argument("--protein-evidence", default="")
    ap.add_argument("--genome-fasta", default="")
    ap.add_argument("--blast-exclude", default="", help='the species left out of the BLAST link, as NCBI writes it, e.g. "Aspergillus fumigatus (taxid:746128)"')
    ap.add_argument("--template", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes, by = read_gff(a.gff)
    pg, occ = read_matrix(a.matrix, a.genome)
    psau = {}
    if a.winners:
        for r in csv.DictReader(open(a.winners), delimiter="\t"):
            psau[r["gene"]] = r.get("psauron", "")

    qual = {}
    if a.quality:
        for r in csv.DictReader(open(a.quality), delimiter="\t"):
            qual[r["pangene"]] = [r.get("strains_valid", ""), r.get("strains_invalid", ""), r.get("strains_not_judged", "")]

    items = []
    for g, x in genes.items():
        if x["new"]:
            p = pg.get(g, "")
            items.append(dict(kind="new", id=g, c=x["c"], s=x["s"], e=x["e"], st=x["st"], pangene=p or "(no pangene)",
                              strains=occ.get(p, 0), model_from=x["model_from"], panel_gene=x["basket"],
                              origin=x["origin"], orf=x["orf"], domain=x["domain"], psauron=psau.get(x["basket"], ""),
                              desc=x["desc"], orfs=qual.get(p)))
    # the window of each item, and the reference genes in it
    for it in items:
        if "c" not in it:
            continue
        # the window: the gene, and the whole gene before and after it (at least FLANK bp each side, at most MAXSIDE)
        L = by.get(it["c"], [])
        left = [x for x in L if x[1] < it["s"]]
        right = [x for x in L if x[0] > it["e"]]
        ws = min(it["s"] - FLANK, max((x for x in left), key=lambda x: x[1])[0]) if left else it["s"] - FLANK
        we = max(it["e"] + FLANK, min((x for x in right), key=lambda x: x[0])[1]) if right else it["e"] + FLANK
        it["win"] = [max(1, ws, it["s"] - MAXSIDE), min(we, it["e"] + MAXSIDE)]
        it["genes"] = [[g, genes[g]["s"], genes[g]["e"], genes[g]["st"], genes[g]["cds"], pg.get(g, ""), genes[g]["desc"]]
                       for (s, e, g) in by.get(it["c"], []) if e >= it["win"][0] and s <= it["win"][1]]

    have_j = bool(a.junctions)
    if have_j:
        wins = defaultdict(list)
        for k, it in enumerate(items):
            if "win" in it:
                wins[it["c"]].append((it["win"][0], it["win"][1], k))
                it["junctions"] = []
        for l in openf(a.junctions):
            f = l.split("\t", 9)
            if len(f) < 9 or f[0] not in wins:
                continue
            s, e = int(f[3]), int(f[4])
            for ws, we, k in wins[f[0]]:
                if e >= ws and s <= we:
                    items[k]["junctions"].append([s, e, f[6], int(float(f[5])), float(att(f[8], "isrpm") or 0),
                                                  att(f[8], "evidence"), att(f[8], "annotatedintron")])
        for it in items:
            if it["kind"] == "new":
                segs = genes[it["id"]]["cds"]
                introns = [(segs[i][1] + 1, segs[i + 1][0] - 1) for i in range(len(segs) - 1)]
                ok = sum(1 for (s, e) in introns if any(abs(j[0] - s) <= 3 and abs(j[1] - e) <= 3 for j in it["junctions"]))
                it["introns"] = [len(introns), ok]

    have_r, baseline = bool(a.region_rnaseq), None
    if have_r:
        reg = {}
        lines = open(a.region_rnaseq).read().splitlines()
        for l in lines:
            if l.startswith("# baseline"):
                baseline = float(l.split("\t")[2])
        for r in csv.DictReader((l for l in lines if not l.startswith("#")), delimiter="\t"):
            reg[("new" if r["kind"] == "new_gene" else "abs", r["id"])] = r
        for it in items:
            r = reg.get(("new" if it["kind"] == "new" else "abs", it["id"]))
            if r:
                it["rnaseq"] = dict(mean=r["all_samples"], pct=r["percentile_own"], best=r["best_sample"],
                                    best_mean=r["best_mean"], breadth=r["best_breadth"], samples=r["samples_with_reads"])
                it["profile"] = [float(v) for v in r["profile"].split(",") if v]
    if a.rnaseq_calls:
        hdr, calls = None, {}
        for l in open(a.rnaseq_calls):
            if l.startswith("#"):
                continue
            f = l.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            calls[f[0]] = dict(zip(hdr, f))
        for it in items:
            if it["kind"] == "new" and it["id"] in calls:
                it.setdefault("rnaseq", {})["call"] = calls[it["id"]].get("call", "")

    if a.genome_fasta:
        from protein_evidence import fasta, protein   # same folder
        seqs = fasta(a.genome_fasta, {it["c"] for it in items if it["kind"] == "new"})
        for it in items:
            if it["kind"] == "new" and it["c"] in seqs and genes[it["id"]]["cds"]:
                it["aa"] = protein(seqs[it["c"]], genes[it["id"]]["cds"], it["st"])
        print(f"BLAST links: {sum(1 for it in items if 'aa' in it)} proteins from {a.genome_fasta}")
    have_p = bool(a.protein_evidence)
    if have_p:
        pe = {r["gene"]: r for r in csv.DictReader(open(a.protein_evidence), delimiter="\t")}
        for it in items:
            r = pe.get(it["id"]) if it["kind"] == "new" else None
            if r:
                it["prot"] = [r["result"], r["best_hit"], r["species"], r["identity"], r["query_coverage"], r["evalue"],
                              r["title"]]
    data = dict(explorer=a.explorer, genome=a.genome, have_junctions=have_j, have_rnaseq=have_r or bool(a.rnaseq_calls), baseline=baseline,
                have_protein=have_p, blast_exclude=a.blast_exclude, items=items)
    open(a.out, "w").write(open(a.template).read().replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    n = defaultdict(int)
    for it in items:
        n[it["kind"]] += 1
    print(dict(n), "-> " + a.out)


if __name__ == "__main__":
    main()
