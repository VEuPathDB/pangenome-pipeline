#!/usr/bin/env python3
"""
Baskets: genes at the same place, on the same strand, go into the same basket and compete; the reference genes always
win; other models win when they do not clash with a winner; a basket can have several winners. Each winner is one gene
of the panel.

1. Baskets. Nodes are native genes. a (genome A) and b (genome B) are linked when a's Liftoff transfer onto B (or b's
   onto A) overlaps b (or a) on the same strand by >= --link-bp of CDS span; two native genes of the same genome are
   linked when they overlap that much. Genes on opposite strands are never linked. No reciprocity: pieces, fusions and
   duplicates join and compete.
2. Place check. Each member is placed on the reference genome: its own CDS span (reference genes) or its primary
   Liftoff transfer onto the reference (coverage >= --min-ref-cov). Members whose places overlap (same strand,
   >= --link-bp) stay together; the group in the most genomes keeps the basket, any other group (a copy at another
   place, e.g. a paralogue) forms a basket of its own. A member with no place on the reference is checked the same
   way on the next genome of --genomes where it has a place (its own gene there, or its transfer with coverage
   >= --min-ref-cov) and where the members that keep the basket have one too: it stays when its place overlaps
   theirs, else it forms a basket of its own. Only when no genome places it is it judged by its neighbours (baskets
   of the nearest --k genes on each side): it leaves only when it has >= --min-shared neighbour baskets and shares
   fewer than --min-shared of them with the members that keep the basket.
   Copies (checked first, in each basket): two members of the same genome at two places of that genome (not
   overlapping, and not a tandem array: another contig, or more than --array-k genes apart) whose proteins align
   (blastp: identity >= --array-id over >= --array-cov of the shorter) are two copies, not one gene: the copy that
   shares fewer neighbour baskets with the other genomes' members leaves and forms a basket of its own
   (copies_split.tsv), and is judged with --psauron-min (the reference annotates the gene). Pieces of a broken gene are different parts of it, their proteins do not align, and they stay.
   A reference genome's own genes are never split.
3. ORFs are checked on the genome sequence (start, stop, frame, no internal stop, no N; check_orf.py), for every
   member's own model and every transfer. Selenoproteins (description) count as valid.
4. Competition, per basket. model(s, t) = member s's own model in its genome, its primary transfer elsewhere.
   Every reference member wins. The others are ranked by:
     0. no N (assembly gap) in its own CDS: a model with an N is never preferred;
     a. protein domain (has_domain of the member, from --domain-dir), unless the member is a fusion (its place on the
        reference overlaps >= 2 reference genes);
     b. PSAURON in-frame score of its own model (only with --psauron-dir or --run-psauron);
     c. valid votes: genomes whose own member agrees with model(s, t) (CDS spans overlap >= --agree of each) and has
        a valid ORF, +1 for its own genome when its own model is valid;
     d. valid transfers: genomes where model(s, t) has a valid ORF and coverage >= --min-cov;
     e. CDS span closest to the median of the basket's members.
   A non-reference member can win only when it has a domain (not a fusion), or when all of these hold:
     - a valid ORF in >= 1 genome (its own model, or a transfer with coverage >= --min-cov);
     - its own protein is >= --min-aa amino acids;
     - (with PSAURON) its score is >= --psauron-min-in-ref when its transfer onto the reference has coverage
       >= --min-ref-cov (the reference has the sequence but does not annotate a gene there), else >= --psauron-min.
   Each eligible member is accepted unless it clashes with an accepted winner (same strand, overlap >= --clash-bp, in
   at least half of the genomes where both have a model). A basket without a reference gene and without an eligible
   member has no winner: its genes are not written (not_eligible.tsv lists why).
5. A winner is present in genome t when model(winner, t) exists with coverage >= --min-cov (or is its own gene).

6. members.tsv: every member, with the quality check of its own model at its own place (used by the build for the
   members that no model of their basket covers): passes with a domain (not a fusion), else needs a valid ORF (an N
   or a model within --contig-end bp of a contig end cannot be judged and counts as valid), >= --min-aa aa and a
   PSAURON score >= --psauron-min (a basket with a reference gene: the reference annotates the gene) or the winner
   thresholds (a basket without one). Empty keep_check = passes; otherwise the reason.

7. Tandem arrays: a basket where a genome has >= 2 copies of its own (non-overlapping, same strand, at most --array-k
   genes apart, proteins >= --array-id % identical over >= --array-cov of each, blastp) holds a tandem array, whose
   copy number can differ between genomes. Written to tandem_arrays.tsv; decide_models.py and the build keep each
   genome's own copies there instead of one model per winner.

PSAURON scores (one of the two, or neither = no PSAURON rule):
  --psauron-dir DIR   scores computed beforehand (e.g. step 04): every *.csv in DIR, PSAURON output format; IDs are
                      gene IDs or transcript IDs of the native GFFs. A candidate without a score stops the run.
  --run-psauron       the script scores the candidates itself with --psauron-cmd: the members of baskets without a
                      reference gene, and the non-reference members of other baskets that do not clash with the
                      reference gene. Proteins in <out-dir>/psauron_candidates.faa, scores cached in
                      <out-dir>/psauron_candidates.csv (re-used: only genes still without a score are run).
"""
import argparse
import bisect
import csv
import glob
import os
import re
import shlex
import statistics
import subprocess
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_baskets import read_gff   # noqa: E402
from decide_models import read_models, spliced   # noqa: E402
import check_orf   # noqa: E402

DESC = re.compile(r"(?:^|;)description=([^;]*)")
CODE = {a + b + c: aa for (a, b, c), aa in zip([(x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"],
                                               "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


def protein(cds):
    """translation of a CDS without stop codons ('*' removed, as in the step 04 protein files)."""
    s = cds.upper()
    return "".join(CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3)).replace("*", "")


def protein_pairs(seqs, by_genome, threads):
    """blastp of the candidate copies of each genome against each other (one search per genome):
    (x, y) -> (identity %, share of x aligned, share of y aligned), best HSP, for the pairs wanted."""
    import tempfile
    best = {}
    for g, pairs in by_genome.items():
        ids = sorted({m for p in pairs for m in p})
        with tempfile.TemporaryDirectory() as tmp:
            fa = os.path.join(tmp, "p.faa")
            with open(fa, "w") as fh:
                for k in ids:
                    fh.write(f">{k}\n{seqs[k]}\n")
            subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "prot", "-out", os.path.join(tmp, "db")], check=True,
                           stdout=subprocess.DEVNULL)
            res = subprocess.run(["blastp", "-query", fa, "-db", os.path.join(tmp, "db"), "-evalue", "1e-5",
                                  "-max_target_seqs", "50", "-num_threads", str(threads), "-outfmt",
                                  "6 qseqid sseqid pident qstart qend sstart send qlen slen bitscore"],
                                 check=True, capture_output=True, text=True).stdout
        for line in res.splitlines():
            f = line.split("\t")
            k = (f[0], f[1])
            if k not in pairs and (f[1], f[0]) not in pairs:
                continue
            if k not in best or float(f[9]) > best[k][0]:
                best[k] = (float(f[9]), float(f[2]), (int(f[4]) - int(f[3]) + 1) / int(f[7]),
                           (int(f[6]) - int(f[5]) + 1) / int(f[8]))
    return {k: v[1:] for k, v in best.items()}


def read_psauron(path):
    """PSAURON output (comment lines first, then description,psauron_is_protein,in-frame_score)
    -> {ID: (psauron_is_protein, in-frame score)}."""
    out, header = {}, False
    for row in csv.reader(open(path)):
        if not header:
            header = bool(row) and row[0] == "description"
            continue
        if len(row) >= 3:
            try:
                out[row[0]] = (row[1], float(row[2]))
            except ValueError:
                pass
    return out


def ov(a, b):
    """overlap of two (contig, strand, s, e) on the same strand."""
    if a[0] != b[0] or a[1] != b[1]:
        return 0
    return max(0, min(a[3], b[3]) - max(a[2], b[2]) + 1)


def ln(a):
    return a[3] - a[2] + 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--transfer-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--domain-dir", required=True, help="<genome>_gene_classification.tsv (has_domain of each gene)")
    ap.add_argument("--reference", default="")
    ap.add_argument("--link-bp", type=int, default=20)
    ap.add_argument("--agree", type=float, default=0.8)
    ap.add_argument("--clash-bp", type=int, default=20)
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--min-ref-cov", type=float, default=0.5, help="a transfer onto the reference places a member")
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--min-shared", type=int, default=2)
    ap.add_argument("--min-aa", type=int, default=100, help="shortest protein a non-reference member can win with")
    ap.add_argument("--psauron-dir", default="", help="PSAURON scores computed beforehand (*.csv)")
    ap.add_argument("--run-psauron", action="store_true", help="score the candidates with --psauron-cmd")
    ap.add_argument("--psauron-cmd", default="psauron", help='e.g. "singularity exec GAQA_psauron.sif psauron"')
    ap.add_argument("--psauron-min", type=float, default=0.7, help="minimum score when the reference lacks the sequence")
    ap.add_argument("--psauron-min-in-ref", type=float, default=0.9,
                    help="minimum score when the member's transfer onto the reference works")
    ap.add_argument("--array-k", type=int, default=3,
                    help="tandem arrays: two copies of a genome are at most this many genes apart")
    ap.add_argument("--array-id", type=float, default=70.0, help="tandem arrays: minimum protein identity of two copies")
    ap.add_argument("--array-cov", type=float, default=0.7,
                    help="tandem arrays: minimum share of each protein covered by the alignment of two copies")
    ap.add_argument("--contig-end", type=int, default=1000,
                    help="members.tsv: an own model this close to a contig end counts as a valid ORF")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--check", default="")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    if a.psauron_dir and a.run_psauron:
        ap.error("use either --psauron-dir or --run-psauron, not both")
    use_psauron = bool(a.psauron_dir or a.run_psauron)
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    ref = a.reference or genomes[0]

    native, owner, nmodel, seleno = {}, {}, {}, set()
    tx2gene = {}   # transcript ID -> gene ID (PSAURON files may be keyed by transcript)
    idx = defaultdict(list)
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))
        for line in open(path):
            f = line.split("\t")
            if len(f) > 8 and f[2] in ("mRNA", "transcript"):
                i, p = re.search(r"ID=([^;\n]+)", f[8]), re.search(r"Parent=([^;\n]+)", f[8])
                if i and p:
                    tx2gene[i.group(1)] = p.group(1).split(",")[0]
        genes = read_gff(path)
        for gid, (c, st, s, e, attr) in genes.items():
            native[gid] = (c, st, s, e)
            owner[gid] = g
            idx[(g, c, st)].append((s, e, gid))
            d = DESC.search(attr)
            if d and "selenoprotein" in d.group(1).lower():
                seleno.add(gid)
        nmodel.update(read_models(path, set(genes)))
    for k in idx:
        idx[k].sort()
    starts = {k: [x[0] for x in v] for k, v in idx.items()}
    dom = {}
    for g in genomes:
        for r in csv.DictReader(open(os.path.join(a.domain_dir, f"{g}_gene_classification.tsv")), delimiter="\t"):
            dom[r["gene_id"]] = r["has_domain"] == "True"

    def natives_over(g, m):
        """native genes of genome g overlapping model m on the same strand by >= --link-bp."""
        c, st, s, e = m[:4]
        L, S = idx.get((g, c, st), []), starts.get((g, c, st), [])
        i = bisect.bisect_right(S, e)
        out = []
        for j in range(i - 1, -1, -1):
            bs, be, b = L[j]
            if be < s - 300000:
                break
            if min(e, be) - max(s, bs) + 1 >= a.link_bp:
                out.append(b)
        return out

    parent = {g: g for g in native}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        x, y = find(x), find(y)
        if x != y:
            parent[x] = y

    # 1. links: same-genome overlaps, then transfers
    for k, L in idx.items():
        for i, (s, e, gid) in enumerate(L):
            for s2, e2, g2 in L[i + 1:i + 30]:
                if s2 > e:
                    break
                if min(e, e2) - s2 + 1 >= a.link_bp:
                    union(gid, g2)
    tr, tmodel = {}, {}
    for path in sorted(glob.glob(os.path.join(a.transfer_dir, "*_onto_*"))):
        src, tgt = os.path.basename(path).split("_liftoff")[0].split("_onto_")
        if src not in genomes or tgt not in genomes:
            continue
        want = {k for k in native if owner[k] == src}
        for gid, m in read_models(path, want, primary_only=True).items():
            c, st, segs, cov = m
            tr[(gid, tgt)] = (c, st, segs[0][0], segs[-1][1], None, cov)
            tmodel[(gid, tgt)] = m
            for b in natives_over(tgt, tr[(gid, tgt)]):
                union(gid, b)
    print(f"{len(native)} native genes, {len(tr)} transfers", file=sys.stderr)

    # 3. ORFs on the genome sequence
    own_valid, own_gap, own_end, prot = {}, set(), set(), {}
    for g in genomes:
        genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
        for gid, m in nmodel.items():
            if owner.get(gid) == g:
                if m[2][0][0] <= a.contig_end or m[2][-1][1] > len(genome.get(m[0], "")) - a.contig_end:
                    own_end.add(gid)
                seq = spliced(genome, m)
                own_valid[gid] = gid in seleno or bool(check_orf.check(seq)[-1])
                if "N" in seq.upper():
                    own_gap.add(gid)
                if g != ref:
                    prot[gid] = protein(seq)
        for (gid, t), m in tmodel.items():
            if t == g:
                c, st, s, e, _, cov = tr[(gid, t)]
                tr[(gid, t)] = (c, st, s, e, gid in seleno or bool(check_orf.check(spliced(genome, m))[-1]), cov)
        del genome
        print(f"  ORFs of {g}: done", file=sys.stderr)
    del tmodel

    baskets = defaultdict(list)
    for gid in native:
        baskets[find(gid)].append(gid)

    # 2. place check on the reference, neighbours as a fallback
    basket_of = {m: b for b, ms in baskets.items() for m in ms}
    order = defaultdict(list)
    for gid, (c, st, s, e) in native.items():
        order[(owner[gid], c)].append((s, gid))
    pos = {}
    for k in order:
        order[k].sort()
        for i, (_, gid) in enumerate(order[k]):
            pos[gid] = (k, i)

    def neighbourhood(gid):
        k, i = pos[gid]
        L, own = order[k], basket_of[gid]
        out = set()
        for step in (-1, 1):
            j, got = i + step, []
            while 0 <= j < len(L) and len(got) < a.k:
                b = basket_of[L[j][1]]
                if b != own and b not in got:
                    got.append(b)
                j += step
            out.update(got)
        return out

    def place_on(m, g):
        """member m's place on genome g: its own CDS span there, or its transfer onto g (coverage >= --min-ref-cov)."""
        if owner[m] == g:
            return native[m]
        t = tr.get((m, g))
        return t[:4] if t and t[5] >= a.min_ref_cov else None

    def place(m):
        return place_on(m, ref)

    # copies: two members of one genome at two places of that genome (not overlapping, not a tandem array: another
    # contig or more than --array-k genes apart) whose proteins align (blastp, identity >= --array-id over >= --array-cov
    # of the shorter). Pieces of a broken gene are different parts of it: their proteins do not align.
    pair_by_g, pairs_in = defaultdict(set), defaultdict(list)
    for b, ms in baskets.items():
        by_g = defaultdict(list)
        for m in ms:
            if owner[m] != ref and prot.get(m):
                by_g[owner[m]].append(m)
        for g, gm in by_g.items():
            for i, x in enumerate(gm):
                for y in gm[i + 1:]:
                    (kx, ix), (ky, iy) = pos[x], pos[y]
                    if ov(native[x], native[y]) > 0 or (kx == ky and abs(ix - iy) <= a.array_k):
                        continue
                    pair_by_g[g].add((x, y))
                    pairs_in[b].append((x, y))
    print(f"copies: {sum(len(v) for v in pair_by_g.values())} pairs of members of one genome at two places to compare",
          file=sys.stderr, flush=True)
    copy_sim = protein_pairs(prot, pair_by_g, a.threads)

    def is_copy(x, y):
        sc = copy_sim.get((x, y))
        if sc is None:
            sc = copy_sim.get((y, x))
            sc = sc and (sc[0], sc[2], sc[1])
        if not sc:
            return False
        short = sc[1] if len(prot[x]) <= len(prot[y]) else sc[2]
        return sc[0] >= a.array_id and short >= a.array_cov

    def split_copies(b, ms):
        """the members of a basket that are a second copy (see Copies above): (members that stay, [(copy, kept copy)])."""
        nb_others = {}

        def synteny(x):
            g = owner[x]
            if g not in nb_others:
                nb_others[g] = set()
                for o in ms:
                    if owner[o] != g:
                        nb_others[g] |= neighbourhood(o)
            return len(neighbourhood(x) & nb_others[g])

        left, out_m = [], set()
        for x, y in sorted(pairs_in.get(b, [])):
            if x in out_m or y in out_m or not is_copy(x, y):
                continue
            keep, go = sorted((x, y), key=lambda m: (synteny(m), ln(native[m]), m), reverse=True)
            out_m.add(go)
            left.append((go, keep))
        return [m for m in ms if m not in out_m], left

    final, copy_rows = [], []
    n_split = n_neigh = 0
    n_fallback = Counter()   # unplaced members decided on another genome: stayed / left
    for b, ms in baskets.items():
        if len(ms) < 2:
            final.append(ms)
            continue
        ms, left = split_copies(b, ms)
        for m, other in left:
            final.append([m])
            copy_rows.append([m, owner[m], other])
        if len(ms) < 2:
            final.append(ms)
            continue
        pl = {m: place(m) for m in ms}
        located = [m for m in ms if pl[m]]
        p2 = {m: m for m in located}

        def f2(x):
            while p2[x] != x:
                p2[x] = p2[p2[x]]
                x = p2[x]
            return x

        for i, m1 in enumerate(located):
            for m2 in located[i + 1:]:
                if ov(pl[m1], pl[m2]) >= a.link_bp:
                    p2[f2(m1)] = f2(m2)
        groups = defaultdict(list)
        for m in located:
            groups[f2(m)].append(m)
        groups = sorted(groups.values(), key=lambda g: (-len({owner[x] for x in g}), -len(g)))
        unplaced = [m for m in ms if not pl[m]]
        if not groups:
            final.append(ms)
            continue
        main_g = list(groups[0])
        nb_main = None
        for m in unplaced:
            # the next genome of the list where m and the members that keep the basket both have a place
            stays = None
            for g in genomes:
                if g == ref:
                    continue
                pm = place_on(m, g)
                if not pm:
                    continue
                theirs = [p for p in (place_on(x, g) for x in main_g) if p]
                if theirs:
                    stays = any(ov(pm, p) >= a.link_bp for p in theirs)
                    break
            if stays is None:   # no genome places it: its neighbours decide
                if nb_main is None:
                    nb_main = set()
                    for x in main_g:
                        nb_main |= neighbourhood(x)
                nb = neighbourhood(m)
                stays = len(nb) < a.min_shared or len(nb & nb_main) >= a.min_shared
                if not stays:
                    n_neigh += 1
            else:
                n_fallback["stayed" if stays else "left"] += 1
            if stays:
                main_g.append(m)
            else:
                final.append([m])
        final.append(main_g)
        for gr in groups[1:]:
            final.append(gr)
            n_split += 1

    # 4. competition
    def model(s, t):
        if owner[s] == t:
            return native[s] + (own_valid.get(s, False), 1.0)
        return tr.get((s, t))

    ref_spans = idx

    def fusion(s):
        p = place(s)
        if not p or owner[s] == ref:
            return False
        L, S = ref_spans.get((ref, p[0], p[1]), []), starts.get((ref, p[0], p[1]), [])
        i = bisect.bisect_right(S, p[3])
        hit = 0
        for j in range(i - 1, -1, -1):
            bs, be, _ = L[j]
            if be < p[2] - 300000:
                break
            if min(p[3], be) - max(p[2], bs) + 1 >= a.link_bp:
                hit += 1
        return hit >= 2

    def clashes(s, w):
        both = hit = 0
        for t in genomes:
            m1, m2 = model(s, t), model(w, t)
            if not m1 or not m2:
                continue
            both += 1
            if ov(m1, m2) >= a.clash_bp:
                hit += 1
        return bool(both) and hit * 2 >= both

    # 4a. ranking information of every member, and the members that could win without a reference gene (candidates)
    per_basket, cands = [], set()
    for ms in final:
        by_g = defaultdict(list)
        for m in ms:
            by_g[owner[m]].append(m)
        med = statistics.median(ln(native[m]) for m in ms)
        info = {}
        for s in ms:
            votes, val = (1 if own_valid.get(s) else 0), 0
            for t in genomes:
                md = model(s, t)
                if not md:
                    continue
                if md[4] and md[5] >= a.min_cov:
                    val += 1
                if t != owner[s]:
                    for m in by_g.get(t, []):
                        o = ov(md, native[m])
                        if o >= a.agree * ln(md) and o >= a.agree * ln(native[m]) and own_valid.get(m):
                            votes += 1
                            break
            dok = dom.get(s, False) and not fusion(s)
            info[s] = [int(s not in own_gap), int(dok), 0.0, votes, val, -abs(ln(native[s]) - med)]
        refs = [m for m in ms if owner[m] == ref]
        for s in ms:
            if owner[s] != ref and not any(clashes(s, w) for w in refs):
                cands.add(s)
        per_basket.append((ms, by_g, info, refs))

    # 4b. PSAURON scores of the candidates
    score = {}   # gene -> (psauron_is_protein, in-frame score)
    if a.psauron_dir:
        for f in sorted(glob.glob(os.path.join(a.psauron_dir, "*.csv"))):
            for k, v in read_psauron(f).items():
                gk = k if k in native else tx2gene.get(k, tx2gene.get(re.sub(r"-p\d+$", "", k), k))
                score[gk] = v
        missing = sorted(c for c in cands if c not in score)
        if missing:
            sys.exit(f"ERROR: {len(missing)} candidates have no PSAURON score in {a.psauron_dir} "
                     f"(e.g. {', '.join(missing[:5])}); score them, or run with --run-psauron instead")
    elif a.run_psauron:
        cache = os.path.join(a.out_dir, "psauron_candidates.csv")
        if os.path.exists(cache):
            score.update(read_psauron(cache))
        for c in cands:
            if not prot.get(c):
                score.setdefault(c, ("False", 0.0))   # no protein (empty model): cannot be scored
        todo = sorted(c for c in cands if c not in score)
        print(f"PSAURON: {len(cands)} candidates, {len(cands) - len(todo)} already scored, {len(todo)} to run",
              file=sys.stderr)
        if todo:
            faa = os.path.join(a.out_dir, "psauron_candidates.faa")
            with open(faa, "w") as fh:
                fh.writelines(f">{c}\n{prot[c]}\n" for c in todo)
            new_csv = os.path.join(a.out_dir, "psauron_candidates.new.csv")
            # real absolute paths, run from the real output folder: a container cannot follow a symlinked working
            # directory (e.g. ~/hc-storage -> /mnt/hc-storage/...) and would not find relative paths
            cmd = shlex.split(a.psauron_cmd) + ["-i", os.path.realpath(faa), "-o", os.path.realpath(new_csv), "-p"]
            print("running: " + " ".join(cmd), file=sys.stderr)
            subprocess.run(cmd, check=True, cwd=os.path.realpath(a.out_dir))
            new = read_psauron(new_csv)
            score.update(new)
            os.remove(new_csv)
            with open(cache, "w") as fh:
                fh.write("description,psauron_is_protein,in-frame_score\n")
                fh.writelines(f"{k},{v[0]},{v[1]}\n" for k, v in sorted(score.items()))
            missing = [c for c in todo if c not in new]
            if missing:
                sys.exit(f"ERROR: PSAURON returned no score for {len(missing)} candidates (e.g. {', '.join(missing[:5])})")
    for _, _, info, _ in per_basket:
        for s in info:
            if s in score:
                info[s][2] = score[s][1]

    split_copy = {r[0] for r in copy_rows}

    def not_eligible(s, info):
        """why a non-reference candidate cannot win ("" = it can)."""
        if info[s][1]:
            return ""
        if info[s][4] == 0:
            return "no valid ORF in any genome"
        if len(prot.get(s, "")) < a.min_aa:
            return f"protein shorter than {a.min_aa} aa"
        if use_psauron:
            # a split copy: the reference annotates the gene (the other copy), so --psauron-min
            in_ref = place(s) and s not in split_copy
            thr = a.psauron_min_in_ref if in_ref else a.psauron_min
            if info[s][2] < thr:
                return f"PSAURON {info[s][2]:.3f} < {thr}" + (" (sequence in the reference)" if in_ref else "")
        return ""

    def own_model_check(s, info, has_ref):
        """members.tsv: why member s's own model cannot stand at its own place ("" = it can)."""
        if owner[s] == ref or info[s][1]:
            return ""
        if not (own_valid.get(s) or s in own_gap or s in own_end):
            return "invalid ORF"
        if len(prot.get(s, "")) < a.min_aa:
            return f"protein shorter than {a.min_aa} aa"
        if use_psauron and s in score:
            thr = a.psauron_min if has_ref or s in split_copy else (a.psauron_min_in_ref if place(s) else a.psauron_min)
            if score[s][1] < thr:
                return f"PSAURON {score[s][1]:.3f} < {thr}"
        return ""

    # 4c. competition
    stats = Counter()
    winners_all, rows, wrows, nerows, array_input, mrows = [], [], [], [], [], []
    present = Counter()
    checks = {x for x in a.check.split(",") if x}
    for ms, by_g, info, refs in per_basket:
        mset = set(ms)
        others = sorted((m for m in ms if owner[m] != ref), key=lambda m: (info[m], m), reverse=True)
        win = list(refs)
        for s in others:
            if s in cands:
                why = not_eligible(s, info)
                if why:
                    nerows.append([s, owner[s], ";".join(refs) or "-", why, info[s][1], info[s][4],
                                   len(prot.get(s, "")), f"{info[s][2]:.3f}" if s in score else "",
                                   "yes" if place(s) else "no"])
                    stats["members not eligible: " + ("low PSAURON" if why.startswith("PSAURON") else why)] += 1
                    continue
            if not any(clashes(s, w) for w in win):
                win.append(s)
        if not win:
            stats["baskets without a winner (no eligible member): genes not written"] += 1
        stats["baskets"] += 1
        stats[f"winners per basket: {min(len(win), 4) if len(win) < 4 else '4+'}"] += 1
        stats["panel genes (winners)"] += len(win)
        stats["panel genes with a reference gene"] += sum(1 for w in win if owner[w] == ref)
        stats["panel genes without a reference gene"] += sum(1 for w in win if owner[w] != ref)
        for w in win:
            n = 0
            for t in genomes:
                md = model(w, t)
                if md and (owner[w] == t or md[5] >= a.min_cov):
                    present[t] += 1
                    n += 1
            # opposite-strand overlap with a reference gene, in the reference genome
            anti = ""
            mr = model(w, ref)
            if mr and owner[w] != ref:
                c, st, s0, e0 = mr[:4]
                other = "-" if st == "+" else "+"
                for s1, e1, gid in idx.get((ref, c, other), []):
                    if min(e0, e1) - max(s0, s1) + 1 >= 20:
                        anti = gid
                        break
            wrows.append([w, owner[w], int(owner[w] == ref), info[w][3], info[w][4], len(by_g), n, ln(native[w]), anti,
                          info[w][1], int(bool(own_valid.get(w))), len(prot[w]) if w in prot else "",
                          f"{score[w][1]:.3f}" if w in score else ""])
            stats[f"panel genes in {'15' if n == len(genomes) else '2-' + str(len(genomes) - 1) if n > 1 else '1'} genomes"] += 1
        winners_all.extend(win)
        rows.append([";".join(sorted(win)), len(ms), ",".join(sorted(ms))])
        for s in ms:
            mrows.append([s, owner[s], ";".join(sorted(win)) or "-", int(owner[s] == ref), len(prot.get(s, "")),
                          "valid" if own_valid.get(s) else "gap" if s in own_gap else "end" if s in own_end else "invalid",
                          info[s][1], f"{score[s][1]:.3f}" if s in score else "", own_model_check(s, info, bool(refs))])
        if win:
            array_input.append((";".join(sorted(win)), ms))
        hit_checks = checks & mset
        if hit_checks:
            print(f"\n[check] basket with {sorted(hit_checks)}: {len(ms)} members from {len(by_g)} genomes")
            print(f"        winners: {sorted(win)}")
            for m in sorted(ms, key=lambda m: info[m], reverse=True)[:8]:
                print(f"          {m:18} no N {info[m][0]}  domain {info[m][1]}  PSAURON {info[m][2]:.2f}"
                      f"  valid votes {info[m][3]:2}  valid transfers {info[m][4]:2}  span {ln(native[m])}")
    ref_genes = [g for g in native if owner[g] == ref]
    lost_ref = [g for g in ref_genes if g not in set(winners_all)]

    # 5. tandem arrays: a basket where a genome has >= 2 copies of its own (non-overlapping, same strand, at most
    # --array-k genes apart, proteins >= --array-id identical over >= --array-cov of each). There the copy number can
    # differ between genomes, so the next steps keep each genome's own copies instead of one model per winner.
    pairs = []
    for key, ms in array_input:
        by_genome = defaultdict(list)
        for m in ms:
            by_genome[owner[m]].append(m)
        for g, gm in by_genome.items():
            for i, x in enumerate(gm):
                for y in gm[i + 1:]:
                    (kx, ix), (ky, iy) = pos[x], pos[y]
                    if kx != ky or native[x][1] != native[y][1] or ov(native[x], native[y]) > 0:
                        continue
                    lx, ly = len(prot.get(x, "")), len(prot.get(y, ""))
                    # each copy must cover >= --array-cov of the other: lengths within that ratio (pieces are left out)
                    if abs(ix - iy) <= a.array_k and lx and ly and min(lx, ly) >= a.array_cov * max(lx, ly):
                        pairs.append((key, g, x, y))
    by_genome = defaultdict(set)
    for key, g, x, y in pairs:
        by_genome[g].add((x, y))
    print(f"tandem arrays: {len(pairs)} neighbouring pairs to compare, "
          f"{len({m for v in by_genome.values() for p in v for m in p})} proteins", file=sys.stderr, flush=True)
    sim = protein_pairs(prot, by_genome, a.threads)
    copies = defaultdict(lambda: defaultdict(set))   # basket key -> genome -> copies
    for key, g, x, y in pairs:
        sc = sim.get((x, y)) or sim.get((y, x))
        if sc and sc[0] >= a.array_id and sc[1] >= a.array_cov and sc[2] >= a.array_cov:
            copies[key][g].update((x, y))
    with open(os.path.join(a.out_dir, "tandem_arrays.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["winners", "copies_per_genome", "copies"])
        for key in sorted(copies):
            w.writerow([key, ";".join(f"{g}:{len(c)}" for g, c in sorted(copies[key].items())),
                        ",".join(sorted(m for c in copies[key].values() for m in c))])
    print(f"\ntandem arrays: {len(copies)} baskets with >= 2 copies in a genome ({len(pairs)} neighbouring pairs "
          f"compared)")
    with open(os.path.join(a.out_dir, "baskets_v2.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["winners", "n_members", "members"])
        w.writerows(rows)
    with open(os.path.join(a.out_dir, "winners.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "genome", "is_reference", "support", "valid_orf_genomes", "annotated_in", "present_in",
                    "cds_span", "antisense_to_reference_gene", "has_domain_not_fusion", "own_valid_orf",
                    "protein_aa", "psauron"])
        w.writerows(wrows)
    with open(os.path.join(a.out_dir, "members.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "genome", "basket_winners", "is_reference", "protein_aa", "own_orf",
                    "has_domain_not_fusion", "psauron", "keep_check"])
        w.writerows(mrows)
    with open(os.path.join(a.out_dir, "copies_split.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["copy", "genome", "kept_copy"])
        w.writerows(copy_rows)
    with open(os.path.join(a.out_dir, "not_eligible.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "genome", "reference_genes_in_basket", "reason", "has_domain_not_fusion",
                    "valid_orf_genomes", "protein_aa", "psauron", "transfer_onto_reference"])
        w.writerows(nerows)
    print(f"\nbaskets split by place on the reference: {n_split}; members without a place on the reference decided on "
          f"another genome: {n_fallback['stayed']} stayed, {n_fallback['left']} left; members split off by neighbours "
          f"(no place on any genome): {n_neigh}; copies split off into baskets of their own: {len(copy_rows)}")
    for k in sorted(stats):
        print(f"  {stats[k]:7}  {k}")
    print(f"  reference genes: {len(ref_genes)}, lost: {len(lost_ref)}")
    print("genes per genome (winners present):", ", ".join(f"{t} {present[t]}" for t in genomes))


if __name__ == "__main__":
    main()
