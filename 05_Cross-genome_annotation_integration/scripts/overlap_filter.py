#!/usr/bin/env python3
"""
Last step of the basket annotation: the overlap rule, on the GFFs as they are after the build and miniprot, and the
final GFFs.

Only OPPOSITE-strand overlaps count (same-strand ones are settled before: the clash in baskets_v2.py and "one place,
one gene" in build_basket_annotation.py). Two panel genes overlap in a genome when their CDS spans (first to last CDS
base, introns included) overlap by >= --min-bp on opposite strands. The pair rules decide who loses in that genome:
  both reference genes                     -> neither
  reference gene vs gene without a domain  -> the other one
  reference gene vs gene with a domain     -> neither
  only one with a domain                   -> the one without
  both with a domain                       -> neither
  neither with a domain                    -> lower PSAURON score (when both panel genes have one, from --winners), then
                                              fewer genomes with a valid ORF, then fewer genomes, then shorter CDS
(reference gene = no not_in_reference_annotation=True; domain = has_domain=True; valid ORF = orf_status valid /
repaired_miniprot / selenoprotein, or no orf_invalid=True when miniprot did not run.)
A panel gene is removed from every genome when it loses in >= min(--min-strains, genomes where it is present).

When a removed gene had taken the place of another model in a genome ("one place, one gene" in the build), that model
is put back from <build-dir>/<genome>_stacked_dropped.gff (restored_after_removal_of=<gene>), unless its own panel
gene is removed too or it overlaps a remaining gene on the same strand by >= 50% of the smaller.

A pseudogene's CDS line (its span, added by prepare_pseudogenes.py) is not written: it ends as gene + transcript + exon.

Writes in --out-dir: <genome>_basket_annotation.gff (final), removed_panel_genes.tsv (the build's list + overlap),
overlap_losses.tsv (every panel gene that loses at least once), restored_genes.tsv.
"""
import argparse
import csv
import os
import re
import statistics
from collections import Counter, defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
BASKET_RE = re.compile(r"(?:^|;)basket=([^;]+)")
STATUS_RE = re.compile(r"(?:^|;)orf_status=([^;]+)")
VALID_STATUS = {"valid", "repaired_miniprot", "selenoprotein"}


def read_blocks(path):
    """the genes of a GFF: list of dicts with the gene's lines, panel gene, contig, strand, CDS span and length."""
    genes, of, order = {}, {}, []
    if not os.path.exists(path):
        return []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if f[2] == "gene" and i:
            g = {"id": i.group(1), "lines": [line], "attrs": f[8], "contig": f[0], "strand": f[6],
                 "start": int(f[3]), "cds": []}
            b = BASKET_RE.search(f[8])
            g["panel"] = b.group(1) if b else i.group(1)
            genes[g["id"]] = g
            of[g["id"]] = g["id"]
            order.append(g["id"])
        elif p:
            gid = of.get(p.group(1).split(",")[0])
            if gid is None:
                continue
            if i:
                of[i.group(1)] = gid
            genes[gid]["lines"].append(line)
            if f[2] == "CDS":
                genes[gid]["cds"].append((int(f[3]), int(f[4])))
    out = []
    for gid in order:
        g = genes[gid]
        if not g["cds"]:
            continue
        g["s"], g["e"] = min(a for a, _ in g["cds"]), max(b for _, b in g["cds"])
        g["len"] = sum(b - a + 1 for a, b in g["cds"])
        st = STATUS_RE.search(g["attrs"])
        g["valid"] = (st.group(1) in VALID_STATUS) if st else "orf_invalid=True" not in g["attrs"]
        out.append(g)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff-dir", required=True, help="GFFs after the build (and miniprot): <genome>_basket_annotation.gff")
    ap.add_argument("--build-dir", required=True, help="build_basket_annotation.py output (stacked_dropped, removed lists)")
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--min-bp", type=int, default=20)
    ap.add_argument("--min-strains", type=int, default=2)
    ap.add_argument("--winners", default="", help="baskets_v2.py's winners.tsv: the PSAURON score of each panel gene")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]

    by_genome = {g: read_blocks(os.path.join(a.gff_dir, f"{g}_basket_annotation.gff")) for g in genomes}

    # per panel gene: reference, domain, genomes present, genomes with a valid ORF, median CDS length
    P = defaultdict(lambda: {"ref": True, "dom": False, "present": 0, "valid": 0, "lens": []})
    for g, genes in by_genome.items():
        for x in genes:
            p = P[x["panel"]]
            p["ref"] = p["ref"] and "not_in_reference_annotation=True" not in x["attrs"]
            p["dom"] = p["dom"] or "has_domain=True" in x["attrs"]
            p["present"] += 1
            p["valid"] += int(x["valid"])
            p["lens"].append(x["len"])
    for p in P.values():
        p["cds"] = statistics.median(p["lens"])
    psauron = {}
    if a.winners:
        for r in csv.DictReader(open(a.winners), delimiter="\t"):
            if r.get("psauron"):
                psauron[r["gene"]] = float(r["psauron"])
    print(f"PSAURON scores: {len(psauron)} ({sum(1 for k in P if k in psauron)} panel genes)")

    def loser(p, q):
        A, B = P[p], P[q]
        if A["ref"] and B["ref"]:
            return None
        if A["ref"] or B["ref"]:
            o = q if A["ref"] else p
            return None if P[o]["dom"] else o
        if A["dom"] != B["dom"]:
            return q if A["dom"] else p
        if A["dom"]:
            return None
        if p in psauron and q in psauron and psauron[p] != psauron[q]:
            return p if psauron[p] < psauron[q] else q
        for k in ("valid", "present", "cds"):
            if A[k] != B[k]:
                return p if A[k] < B[k] else q
        return None

    # losses per genome (a gene counts once per genome)
    lost_in, against = Counter(), defaultdict(set)
    n_pairs = 0
    for g, genes in by_genome.items():
        losers = set()
        by_contig = defaultdict(list)
        for x in genes:
            by_contig[x["contig"]].append(x)
        for xs in by_contig.values():
            xs.sort(key=lambda x: x["s"])
            for i, x in enumerate(xs):
                for y in xs[i + 1:]:
                    if y["s"] > x["e"]:
                        break
                    if x["strand"] == y["strand"] or x["panel"] == y["panel"]:
                        continue
                    if min(x["e"], y["e"]) - max(x["s"], y["s"]) + 1 < a.min_bp:
                        continue
                    n_pairs += 1
                    lo = loser(x["panel"], y["panel"])
                    if lo:
                        losers.add(lo)
                        against[lo].add(y["panel"] if lo == x["panel"] else x["panel"])
        for lo in losers:
            lost_in[lo] += 1
    removed = {p for p, n in lost_in.items() if n >= min(a.min_strains, P[p]["present"])}

    # put back the models a removed gene had pushed out ("one place, one gene" of the build)
    restored_rows, extra = [], defaultdict(list)
    stacked = os.path.join(a.build_dir, "stacked_dropped.tsv")
    todo = defaultdict(dict)   # genome -> dropped panel gene -> the removed gene that had its place
    if os.path.exists(stacked):
        for r in csv.DictReader(open(stacked), delimiter="\t"):
            if r["kept_panel_gene"] in removed and r["panel_gene_dropped"] not in removed:
                todo[r["genome"]][r["panel_gene_dropped"]] = r["kept_panel_gene"]
    for g, want in todo.items():
        present = {x["panel"] for x in by_genome[g]}
        kept = [x for x in by_genome[g] if x["panel"] not in removed]
        for x in read_blocks(os.path.join(a.build_dir, f"{g}_stacked_dropped.gff")):
            if x["panel"] not in want or x["panel"] in present:
                continue
            clash = any(y["contig"] == x["contig"] and y["strand"] == x["strand"] and
                        min(x["e"], y["e"]) - max(x["s"], y["s"]) + 1 >= 0.5 * min(x["e"] - x["s"] + 1, y["e"] - y["s"] + 1)
                        for y in kept + extra[g])
            if clash:
                continue
            first = x["lines"][0].rstrip("\n") + f";restored_after_removal_of={want[x['panel']]}\n"
            x["lines"] = [first] + x["lines"][1:]
            extra[g].append(x)
            present.add(x["panel"])
            restored_rows.append([g, x["panel"], x["id"], want[x["panel"]]])

    # final GFFs
    for g in genomes:
        genes = [x for x in by_genome[g] if x["panel"] not in removed] + extra[g]
        corder = {}
        for x in by_genome[g] + extra[g]:
            corder.setdefault(x["contig"], len(corder))
        with open(os.path.join(a.out_dir, f"{g}_basket_annotation.gff"), "w") as fh:
            fh.write("##gff-version 3\n")
            for x in sorted(genes, key=lambda x: (corder[x["contig"]], x["start"])):
                # a pseudogene's CDS was only its span (prepare_pseudogenes.py): written as gene + transcript + exon
                fh.writelines(l for l in x["lines"] if "pseudogene_span=True" not in l)

    rows = []
    prev = os.path.join(a.build_dir, "removed_panel_genes.tsv")
    if os.path.exists(prev):
        rows = [[r["panel_gene"], r["reason"]] for r in csv.DictReader(open(prev), delimiter="\t")]
    rows += [[p, "overlap (opposite strand)"] for p in sorted(removed)]
    with open(os.path.join(a.out_dir, "removed_panel_genes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["panel_gene", "reason"])
        w.writerows(rows)
    with open(os.path.join(a.out_dir, "overlap_losses.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["panel_gene", "is_reference", "has_domain", "present_in", "valid_orf_genomes", "lost_in",
                    "removed", "lost_against"])
        for p in sorted(lost_in, key=lambda p: -lost_in[p]):
            w.writerow([p, int(P[p]["ref"]), int(P[p]["dom"]), P[p]["present"], P[p]["valid"], lost_in[p],
                        int(p in removed), ",".join(sorted(against[p]))])
    with open(os.path.join(a.out_dir, "restored_genes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "panel_gene", "gene_id", "restored_after_removal_of"])
        w.writerows(restored_rows)

    print(f"opposite-strand overlapping pairs (all genomes): {n_pairs}; panel genes that lose at least once: "
          f"{len(lost_in)}; removed from every genome (lose in >= min({a.min_strains}, present)): {len(removed)}; "
          f"models put back: {len(restored_rows)}")
    for g in genomes:
        print(f"  {g:10} {len(by_genome[g]):6} -> {len([x for x in by_genome[g] if x['panel'] not in removed]) + len(extra[g]):6} genes")


if __name__ == "__main__":
    main()
