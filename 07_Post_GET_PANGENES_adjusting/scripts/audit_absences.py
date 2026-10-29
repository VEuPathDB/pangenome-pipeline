#!/usr/bin/env python3
"""
Look at every absence in the pangenome - a pangene with no gene in a genome - and say why the gene is missing there.

For each pangene P and each genome g where P has no gene:

1. Expected place, from the neighbours. In every genome h that has P (with one gene), walking away from P's gene on each
   side, the first gene whose pangene also has a gene in g is the anchor on that side (L, R). In g, the genes of L and R
   are looked up; when they sit on the same contig, at most --max-gap genes apart, the stretch between them is the
   expected place. The anchor pair found in most genomes is used. When L and R are the same pangene (P's gene lies
   inside a gene of L, e.g. in an intron), the expected place is that gene of L in g.
     6 contig_edge        an anchor is at the end of a contig (in g, or in the genomes with P: no anchor on one side)
     7 neighbours_split   both anchors are found in g, but on different contigs or far apart (break or rearrangement)
     8 region_absent      no anchor is found in g
2. What sits in the expected place, in g's final annotation. Each gene there is compared with the genes of P in the other
   genomes by blastn on DNA, twice (--align both): the CDS span (first CDS start to last CDS end, introns included,
   UTRs excluded), which still aligns when two models of the same gene differ in their exons, and the spliced CDS,
   which still aligns when the introns differ in length between genomes. DNA, not protein, so a model whose reading
   frame is broken still aligns. All the HSPs of a pair are taken together. It is the same gene when, in either
   alignment, the identity is >= --min-id and the HSPs cover >= --min-cov of the shorter of the two sequences.
     1a same_gene_movable        the same gene, in another pangene Q that has no genome in common with P: can be moved
     1b same_gene_shared         the same gene, in a pangene Q that also has genes in genomes that have P (duplication?)
     1c same_gene_no_pangene     the same gene, left out of every pangene (e.g. CDS < 100 bp)
3. When no gene in the place matches: the Liftoff transfers of P's genes onto g, inside the place.
     2a other_gene_no_transfer   other genes there, and no transfer of P
     2b other_gene_over_transfer other genes there, over a transfer of P (the build kept another model in its place)
     3  transfer_not_used        no gene there, but a transfer of P
     4  empty_place_shorter      no gene, no transfer, and the place is much shorter than in the genomes with P
                                 (< --short-frac of their median length): likely a deletion
     5  empty_place              no gene, no transfer, place of normal length

Outputs (in --out-dir):
  absences.tsv   one row per (pangene, genome) absence: code, category, place, anchors, the gene matching P (with its
                 pangene and alignment), the genes in the place, and the transfer of P found there (source genome, ID,
                 location, coverage, identity, valid ORF)
  summary.tsv    absences per category and genome

Needs blastn and makeblastdb on the PATH; otherwise only the standard library.

Usage:
    python3 audit_absences.py --genes-matrix for_get_pangenes_fixed_pangenes/pangene_matrix_genes.tr.tab \\
        --gff-dir tagged_annotation_fixed --genome-dir genomes \\
        --pairwise-dir polished_nonredundant --out-dir absence_audit
"""
import argparse
import bisect
import csv
import os
import re
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

EMPTY = {"", "-", "NA"}
ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
ORIG_RE = re.compile(r"(?:^|;)original_ID=([^;]+)")
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_fasta(path):
    seqs, name, parts = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
    if name:
        seqs[name] = "".join(parts)
    return seqs


def read_gff(path):
    """gene ID -> dict(orig, contig, strand, S, E, segs) (longest transcript's CDS)."""
    genes, t2g, cds = {}, {}, defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
            if i and not p and f[2].endswith("gene"):
                o = ORIG_RE.search(f[8])
                genes[i.group(1)] = dict(orig=o.group(1) if o else i.group(1), contig=f[0], strand=f[6],
                                         S=int(f[3]), E=int(f[4]))
            elif i and p:
                t2g[i.group(1)] = p.group(1).split(",")[0]
            if f[2] == "CDS" and p:
                cds[p.group(1).split(",")[0]].append((int(f[3]), int(f[4])))
    best = {}
    for t, segs in cds.items():
        g = t2g.get(t)
        if g in genes:
            bp = sum(e - s + 1 for s, e in segs)
            if g not in best or bp > best[g][0]:
                best[g] = (bp, sorted(segs))
    for g, d in genes.items():
        d["segs"] = best.get(g, (0, []))[1]
    return genes


def cds_seq(d, genome, spliced=False):
    """DNA of the gene's CDS span (first CDS start to last CDS end), or of the spliced CDS."""
    seq = genome.get(d["contig"])
    if not seq or not d["segs"]:
        return ""
    if spliced:
        s = "".join(seq[a - 1:b] for a, b in d["segs"])
    else:
        s = seq[d["segs"][0][0] - 1:d["segs"][-1][1]]
    if d["strand"] == "-":
        s = s.translate(COMP)[::-1]
    return s.upper()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_final_annotation_tagged.gff")
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--pairwise-dir", required=True)
    ap.add_argument("--pair-suffix", default="_liftoff.gff3_polished")
    ap.add_argument("--max-gap", type=int, default=10, help="most genes between the two neighbours in the genome checked")
    ap.add_argument("--min-id", type=float, default=90.0, help="blastn identity (%%) for the same gene")
    ap.add_argument("--min-cov", type=float, default=0.8, help="share of the shorter sequence the blastn hit must cover")
    ap.add_argument("--align", choices=["both", "span", "spliced"], default="both",
                    help="DNA aligned: the CDS span, the spliced CDS, or both (the same gene when either passes)")
    ap.add_argument("--short-frac", type=float, default=0.5)
    ap.add_argument("--max-ref-proteins", type=int, default=3, help="genes of P (from other genomes) aligned against")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # matrix
    with open(args.genes_matrix, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    cols = [(i, h.split("_for_get_pangenes")[0]) for i, h in enumerate(rows[0]) if i > 0 and h.strip()]
    genomes = [g for _, g in cols]
    cells, pg_of = {}, {}
    for r in rows[1:]:
        if not r or r[0].startswith("chr"):
            continue
        c = {g: r[i].split(",") for i, g in cols if i < len(r) and r[i] not in EMPTY}
        cells[r[0]] = c
        for g, xs in c.items():
            for x in xs:
                pg_of[(g, x)] = r[0]
    print(f"{len(cells)} pangenes, {len(genomes)} genomes", file=sys.stderr)

    # annotations: genes in order on each contig
    G, order, pos = {}, {}, {}
    for g in genomes:
        G[g] = read_gff(os.path.join(args.gff_dir, g + args.gff_suffix))
        byc = defaultdict(list)
        for x, d in G[g].items():
            byc[d["contig"]].append((d["S"], x))
        order[g] = {}
        for c, L in byc.items():
            L.sort()
            order[g][c] = [x for _, x in L]
            for k, x in enumerate(order[g][c]):
                pos[(g, x)] = (c, k)

    def anchor(h, x, step, g, limit=30):
        """walking from x in genome h (step -1 / +1), the first gene in another pangene that also has a gene in g."""
        c, k = pos[(h, x)]
        L, p = order[h][c], pg_of.get((h, x))
        k += step
        n = 0
        while 0 <= k < len(L) and n < limit:
            q = pg_of.get((h, L[k]))
            if q and q != p:
                n += 1
                if g in cells[q]:
                    return q, L[k]
            k += step
        return None, None

    def votes_for(P, g):
        """anchor pairs (L, R) for P in genome g, and the median length of that stretch in the genomes with P."""
        votes, lengths = Counter(), []
        for h, xs in cells[P].items():
            if len(xs) != 1 or (h, xs[0]) not in pos:
                continue
            (l, lx), (r, rx) = anchor(h, xs[0], -1, g), anchor(h, xs[0], +1, g)
            votes[(l, r)] += 1
            if lx and rx:
                a, b = G[h][lx], G[h][rx]
                if a["contig"] == b["contig"]:
                    lengths.append(max(a["E"], b["E"]) - min(a["S"], b["S"]))
        return votes, statistics.median(lengths) if lengths else None

    absences = []
    for P, c in cells.items():
        for g in genomes:
            if g not in c:
                absences.append((P, g))
    print(f"{len(absences)} absences", file=sys.stderr)

    # 1. expected place
    rec = {}
    for P, g in absences:
        votes, med = votes_for(P, g)
        res = dict(category=None, place="", neighbours="", genes_in_place=[], median_len=med)
        found_any, best = False, None
        for (l, r), n in votes.most_common():
            lg = [x for x in cells.get(l, {}).get(g, []) if (g, x) in pos] if l else []
            rg = [x for x in cells.get(r, {}).get(g, []) if (g, x) in pos] if r else []
            if lg or rg:
                found_any = True
            if l and l == r:
                # P's gene lies inside a gene of L (e.g. in an intron): the place is that gene of L in g
                for a in lg:
                    ca, ka = pos[(g, a)]
                    if best is None:
                        best = (ca, ka, ka, l, r, a, a)
                if best:
                    break
                continue
            for a in lg:
                for b in rg:
                    ca, ka = pos[(g, a)]
                    cb, kb = pos[(g, b)]
                    if ca == cb and 0 < abs(ka - kb) <= args.max_gap + 1:
                        lo, hi = sorted((ka, kb))
                        if best is None or hi - lo < best[2] - best[1]:
                            best = (ca, lo, hi, l, r, a, b)
            if best:
                break
        if best and best[5] == best[6]:
            ca, lo, hi, l, r, a, b = best
            A = G[g][a]
            inside = [x for x in order[g][ca] if x != a and G[g][x]["S"] <= A["E"] and G[g][x]["E"] >= A["S"]]
            res.update(place=f"{ca}:{A['S']}-{A['E']}", neighbours=f"{l}|{r}", contig=ca, S=A["S"], E=A["E"],
                       genes_in_place=inside, place_len=A["E"] - A["S"] + 1)
        elif best:
            ca, lo, hi, l, r, a, b = best
            A, B = G[g][a], G[g][b]
            s, e = min(A["E"], B["E"]), max(A["S"], B["S"])
            res.update(place=f"{ca}:{s}-{e}", neighbours=f"{l}|{r}", contig=ca, S=s, E=e,
                       genes_in_place=order[g][ca][lo + 1:hi],
                       place_len=max(A["E"], B["E"]) - min(A["S"], B["S"]))
        else:
            # why no place: contig edge, split, or nothing found
            edge, split = False, False
            top = votes.most_common(1)[0][0] if votes else (None, None)
            if None in top and top != (None, None):
                edge = True   # no anchor on one side in the genomes with P: P's gene is at a contig end there
            for (l, r), n in votes.most_common(3):
                for q in (l, r):
                    for x in (cells.get(q, {}).get(g, []) if q else []):
                        if (g, x) not in pos:
                            continue
                        cc, kk = pos[(g, x)]
                        if kk == 0 or kk == len(order[g][cc]) - 1:
                            edge = True
                        else:
                            split = True
            res["category"] = "contig_edge" if edge else "neighbours_split" if split else (
                "neighbours_split" if found_any else "region_absent")
            l, r = votes.most_common(1)[0][0] if votes else (None, None)
            res["neighbours"] = f"{l}|{r}"
        rec[(P, g)] = res

    # 2. blastn of the CDS of the genes in place against the CDS of P's genes
    genome_seq = {}
    def gseq(g):
        if g not in genome_seq:
            genome_seq.clear()
            genome_seq[g] = read_fasta(os.path.join(args.genome_dir, g + args.genome_suffix))
        return genome_seq[g]
    queries, subjects = {}, {}
    need_q = defaultdict(set)
    for (P, g), res in rec.items():
        for x in res["genes_in_place"]:
            need_q[g].add(x)
    need_s = defaultdict(set)
    for P in {P for (P, g), res in rec.items() if res["genes_in_place"]}:
        ref = sorted(cells[P], key=lambda h: (h != "ME49", h))[:args.max_ref_proteins]
        for h in ref:
            need_s[h].add(cells[P][h][0])
    modes = {"span": [False], "spliced": [True], "both": [False, True]}[args.align]
    for g in genomes:
        if not need_q[g] and not need_s[g]:
            continue
        seq = gseq(g)
        for spl in modes:
            for x in need_q[g]:
                p = cds_seq(G[g][x], seq, spl)
                if p:
                    queries.setdefault(spl, {})[f"{g}|{x}"] = p
            for x in need_s[g]:
                p = cds_seq(G[g][x], seq, spl)
                if p:
                    subjects.setdefault(spl, {})[f"{g}|{x}"] = p
    genome_seq.clear()

    def union(iv):
        tot, cur = 0, None
        for a, b in sorted(iv):
            if cur and a <= cur[1] + 1:
                cur[1] = max(cur[1], b)
            else:
                if cur:
                    tot += cur[1] - cur[0] + 1
                cur = [a, b]
        return tot + (cur[1] - cur[0] + 1 if cur else 0)

    def passes(h):
        return h[0] >= args.min_id and max(h[1], h[2]) >= args.min_cov

    hits = defaultdict(dict)   # query -> subject -> (pident, qcov, scov); the better of the two alignments when both
    for spl in modes:
        Q, Sj = queries.get(spl, {}), subjects.get(spl, {})
        print(f"blastn ({'spliced CDS' if spl else 'CDS span'}): {len(Q)} genes in place against {len(Sj)} genes of "
              f"the pangenes", file=sys.stderr)
        if not Q or not Sj:
            continue
        with tempfile.TemporaryDirectory() as tmp:
            qf, sf = os.path.join(tmp, "q.fna"), os.path.join(tmp, "s.fna")
            with open(qf, "w") as fh:
                for k, v in Q.items():
                    fh.write(f">{k}\n{v}\n")
            with open(sf, "w") as fh:
                for k, v in Sj.items():
                    fh.write(f">{k}\n{v}\n")
            subprocess.run(["makeblastdb", "-in", sf, "-dbtype", "nucl", "-out", os.path.join(tmp, "db")],
                           check=True, stdout=subprocess.DEVNULL)
            out = subprocess.run(["blastn", "-task", "blastn", "-dust", "no", "-query", qf, "-db", os.path.join(tmp, "db"),
                                  "-evalue", "1e-5", "-max_target_seqs", "500", "-num_threads", str(args.threads),
                                  "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send"],
                                 check=True, capture_output=True, text=True).stdout
        # all HSPs of a pair together (an indel splits the alignment): coverage = union of the HSPs on each sequence,
        # identity = mean weighted by HSP length
        hsps = defaultdict(list)
        for line in out.splitlines():
            q, sj, pid, ln, ql, sl, qs, qe, ss, se = line.split("\t")
            hsps[(q, sj)].append((float(pid), int(ln), int(ql), int(sl), int(qs), int(qe),
                                  min(int(ss), int(se)), max(int(ss), int(se))))
        for (q, sj), L in hsps.items():
            ql, sl = L[0][2], L[0][3]
            pid = sum(h[0] * h[1] for h in L) / sum(h[1] for h in L)
            h = (pid, union([(x[4], x[5]) for x in L]) / ql, union([(x[6], x[7]) for x in L]) / sl)
            old = hits[q].get(sj)
            if old is None or (passes(h), max(h[1], h[2]), h[0]) > (passes(old), max(old[1], old[2]), old[0]):
                hits[q][sj] = h

    for (P, g), res in rec.items():
        if res["category"]:
            continue
        ref_ids = {f"{h}|{cells[P][h][0]}" for h in cells[P]}
        best_gene = None
        details = []
        for x in res["genes_in_place"]:
            Q = pg_of.get((g, x), "")
            hs = [(v, s) for s, v in hits.get(f"{g}|{x}", {}).items() if s in ref_ids]
            top = max(hs, key=lambda t: (max(t[0][1], t[0][2]), t[0][0])) if hs else None
            same = bool(top and top[0][0] >= args.min_id and max(top[0][1], top[0][2]) >= args.min_cov)
            details.append(f"{x}({G[g][x]['orig']};{Q or 'no_pangene'};"
                           + (f"{top[0][0]:.0f}%/{top[0][1]:.2f}/{top[0][2]:.2f}" if top else "no_hit") + ")")
            if same and (best_gene is None or max(top[0][1], top[0][2]) > best_gene[2]):
                best_gene = (x, Q, max(top[0][1], top[0][2]), top[0][0])
        res["gene_details"] = ";".join(details)
        if best_gene:
            x, Q = best_gene[0], best_gene[1]
            if not Q:
                res["category"] = "same_gene_no_pangene"
            elif set(cells.get(Q, {})) & set(cells[P]):
                res["category"] = "same_gene_shared"
            else:
                res["category"] = "same_gene_movable"
            res.update(match=x, match_pangene=Q, match_identity=f"{best_gene[3]:.1f}", match_cov=f"{best_gene[2]:.2f}")
        elif res["genes_in_place"]:
            res["category"] = "other_gene"

    # 3. Liftoff, only for the places with no matching gene
    todo = defaultdict(list)
    for (P, g), res in rec.items():
        if res["category"] in (None, "other_gene"):
            todo[g].append((P, res))
    for g, items in todo.items():
        want = defaultdict(set)   # source genome -> original IDs of P's genes
        for P, res in items:
            for h, xs in cells[P].items():
                for x in xs:
                    want[h].add(G[h][x]["orig"] if x in G[h] else x)
        tr = defaultdict(list)     # (source, original ID) -> [(contig, S, E, coverage, identity, valid)]
        for h, ids in want.items():
            path = os.path.join(args.pairwise_dir, f"{h}_onto_{g}{args.pair_suffix}")
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                for line in fh:
                    f = line.split("\t")
                    if len(f) < 9 or not f[2].endswith("gene") or "Parent=" in f[8]:
                        continue
                    i = ID_RE.search(f[8].rstrip("\n"))
                    if i and i.group(1) in ids:
                        a = f[8]
                        cov = re.search(r"coverage=([^;]+)", a)
                        idt = re.search(r"sequence_ID=([^;]+)", a)
                        v = re.search(r"valid_ORFs=([^;\n]+)", a)
                        tr[(h, i.group(1))].append((f[0], int(f[3]), int(f[4]), cov.group(1) if cov else "",
                                                   idt.group(1) if idt else "", (v.group(1) if v else "0") != "0",
                                                   f[6]))
        for P, res in items:
            if "contig" not in res:
                continue
            inplace = []
            for h, xs in cells[P].items():
                for x in xs:
                    o = G[h][x]["orig"] if x in G[h] else x
                    for (c, s, e, cov, idt, v, st) in tr.get((h, o), []):
                        if c == res["contig"] and s <= res["E"] and e >= res["S"]:
                            inplace.append((h, o, c, s, e, cov, idt, v, st))
            if inplace:
                h, o, c, s, e, cov, idt, v, st = max(inplace, key=lambda t: (t[7], float(t[5] or 0) * float(t[6] or 0)))
                res["transfer"] = dict(source=h, id=o, contig=c, start=s, end=e, strand=st, coverage=cov,
                                       identity=idt, valid_orf=v)
                res["category"] = "other_gene_over_transfer" if res["category"] == "other_gene" else "transfer_not_used"
            elif res["category"] == "other_gene":
                res["category"] = "other_gene_no_transfer"
            elif res["category"] is None:
                med = res["median_len"]
                res["category"] = ("empty_place_shorter" if med and res["place_len"] < args.short_frac * med
                                   else "empty_place")

    # outputs
    code = {"same_gene_movable": "1a", "same_gene_shared": "1b", "same_gene_no_pangene": "1c",
            "other_gene_no_transfer": "2a", "other_gene_over_transfer": "2b", "transfer_not_used": "3",
            "empty_place_shorter": "4", "empty_place": "5", "contig_edge": "6", "neighbours_split": "7",
            "region_absent": "8"}
    tcols = ["source", "id", "contig", "start", "end", "strand", "coverage", "identity", "valid_orf"]
    with open(os.path.join(args.out_dir, "absences.tsv"), "w") as fh:
        fh.write("pangene\tgenome\tgenomes_with_pangene\tcode\tcategory\tneighbours\tplace\tplace_length\t"
                 "median_length_elsewhere\tmatch\tmatch_pangene\tmatch_identity\tmatch_cov\tgenes_in_place\t"
                 + "\t".join("transfer_" + c for c in tcols) + "\n")
        for (P, g), res in sorted(rec.items()):
            t = res.get("transfer") or {}
            fh.write("\t".join(map(str, [P, g, len(cells[P]), code[res["category"]], res["category"], res["neighbours"],
                                         res["place"], res.get("place_len", ""), res["median_len"] or "",
                                         res.get("match", ""), res.get("match_pangene", ""),
                                         res.get("match_identity", ""), res.get("match_cov", ""),
                                         res.get("gene_details", "")] + [t.get(c, "") for c in tcols])) + "\n")
    cats = list(code)
    cnt = Counter((res["category"], g) for (P, g), res in rec.items())
    with open(os.path.join(args.out_dir, "summary.tsv"), "w") as fh:
        fh.write("code\tcategory\t" + "\t".join(genomes) + "\ttotal\n")
        for k in cats:
            row = [cnt[(k, g)] for g in genomes]
            fh.write(f"{code[k]}\t{k}\t" + "\t".join(map(str, row)) + f"\t{sum(row)}\n")
    tot = Counter(res["category"] for res in rec.values())
    print("; ".join(f"{code[k]} {k} {tot[k]}" for k in cats) + f" -> {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
