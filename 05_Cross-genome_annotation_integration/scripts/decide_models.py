#!/usr/bin/env python3
"""
Per-genome decision for every gene of the panel (the winners of baskets_v2.py), in two modes at once:
  keep_similar        keep a genome's own model when it is similar to the target gene and has a valid ORF
                      (fewest changes to the existing annotations)
  prefer_panel_model  every genome takes the valid model most similar to the target (the target's own transfer on a tie)

Target of a panel gene = the winner's own model in its own genome (for a reference winner, the reference gene itself, not its
transfers). In genome t, the target model is the winner's own model (t = its genome) or its primary Liftoff transfer.

Models of a panel gene in genome t (the candidates): the target's transfer, t's own member genes of the basket, and the
transfers onto t of the other members (Liftoff coverage >= --min-tcov-rescue). A candidate counts when it lies at the
target's place (over the target's transfer when there is one; otherwise over one of t's own members, or anywhere when t has
none) and is not a fusion (an own gene under 2+ panel genes). Each one is compared with the target CDS (minimap2, the best
collinear alignment of the pair; blastn's best single HSP when minimap2 finds none or a CDS is < --short-cds bp):
coverage of the target CDS, identity, and its length / the target's.
"Valid" here means a valid ORF, or an ORF that cannot be judged: an N in the CDS (assembly gap) or a model within
--contig-end bp of a contig end. Those keep their orf_status later and are not sent to miniprot.

Decisions, per panel gene and genome (prefer_panel_model; the winner's own genome always keeps its own model):
  replace / fill    the most similar valid candidate, covering >= --min-cov of the target CDS and <= --max-len x its
                    length, is the target's transfer ("replace" when t had its own model, "fill" when it had none).
                    Most similar = highest coverage, then identity, then the length closest to the target's; on a tie,
                    the target's transfer
  alternative       ... the most similar valid candidate is another model (t's own gene, or another member's transfer)
  target_invalid    no such candidate, and the target's transfer has Liftoff coverage >= --min-tcov but no valid ORF
                    (an N in the CDS, a contig end or a broken ORF; an N stretch does not align, so the comparison fails):
                    kept and marked (for miniprot)
  rescued           none of the above, but a valid candidate covers >= --min-cov-rescue of the target CDS: the most
                    similar one (marked)
  keep_divergent    none of the above, but t annotated the gene: t keeps its own model
  absent            none of the above (the reason is written to absent.tsv)
keep_similar takes the same decision, except keep_original: t has exactly one own member gene over the target model, with
a valid ORF, covering >= --keep-min-cov of the target CDS, <= --max-len x its length, identity >= --min-id: kept
unchanged.
ORFs are checked here (start, stop, frame, no internal stop), from the genome sequence, with check_orf.py's code.

Tandem arrays (--arrays, baskets_v2.py's tandem_arrays.tsv): their baskets are decided like any other, and besides
the model of each winner a genome keeps its own confirmed copies and gets the transfers of the other genomes' copies
that have a valid ORF and land where it has no gene of the basket (Liftoff coverage >= --min-tcov). Written to
array_models.tsv; the build adds them where they do not overlap the models already written.

Writes decisions.tsv (with the ORF status of the prefer_panel_model choice: valid / gap / end / invalid), absent.tsv,
own_invalid_orf.tsv, array_models.tsv and summary.tsv in --out-dir.
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_orf   # noqa: E402

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
COV_RE = re.compile(r"(?:^|;)coverage=([\d.]+)")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_models(path, want, primary_only=False):
    """gene ID -> (contig, strand, CDS segments of its longest transcript, coverage) for the genes in want."""
    genes, tx, cds = {}, {}, defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
            if i and not p:
                if i.group(1) in want and not (primary_only and "extra_copy_number=" in f[8]
                                               and "extra_copy_number=0" not in f[8]):
                    c = COV_RE.search(f[8])
                    genes[i.group(1)] = (f[0], f[6], float(c.group(1)) if c else 1.0)
            elif i and p and f[2] in ("mRNA", "transcript"):
                tx[i.group(1)] = p.group(1).split(",")[0]
            elif f[2] == "CDS" and p:
                cds[p.group(1).split(",")[0]].append((int(f[3]), int(f[4])))
    out = {}
    for t, segs in cds.items():
        g = tx.get(t, t)
        if g in genes and (g not in out or sum(b - a for a, b in segs) > sum(b - a for a, b in out[g][2])):
            out[g] = (genes[g][0], genes[g][1], sorted(segs), genes[g][2])
    return out


def spliced(genome, m):
    c, st, segs = m[0], m[1], m[2]
    seq = "".join(genome[c][a - 1:b] for a, b in segs) if c in genome else ""
    return seq.translate(COMP)[::-1] if st == "-" else seq


def valid(seq):
    r = check_orf.check(seq)
    return r[-1] if isinstance(r, tuple) else bool(r)


def span(m):
    return (m[0], m[1], m[2][0][0], m[2][-1][1])


def ov(a, b):
    if a[0] != b[0] or a[1] != b[1]:
        return 0
    return max(0, min(a[3], b[3]) - max(a[2], b[2]) + 1)



def orf_status(seq, m, clen, end_bp):
    """valid / gap (an N in the CDS) / end (within end_bp of a contig end) / invalid."""
    if valid(seq):
        return "valid"
    if "N" in seq.upper():
        return "gap"
    c = clen.get(m[0], 0)
    if m[2][0][0] <= end_bp or m[2][-1][1] > c - end_bp:
        return "end"
    return "invalid"


VALID_LIKE = {"valid", "gap", "end"}


def cds_sim(qseqs, sseqs, pairs, threads, minimap2, short_cds):
    """(query, subject) -> (share of the subject CDS aligned (its bases in the alignment, gaps not counted), query length /
    subject length, identity %).
    minimap2: the best alignment of each pair (one collinear chain, never several pieces added up); blastn's best single
    HSP for the pairs minimap2 does not align and for CDS shorter than short_cds."""
    sim = {}
    pairs = {(q, s) for q, s in pairs if qseqs.get(q) and sseqs.get(s)}
    if not pairs:
        return sim
    with tempfile.TemporaryDirectory() as tmp:
        qf, sf = os.path.join(tmp, "q.fna"), os.path.join(tmp, "s.fna")
        with open(qf, "w") as fq:
            for k in {p[0] for p in pairs}:
                fq.write(f">{k}\n{qseqs[k]}\n")
        with open(sf, "w") as fs:
            for k in {p[1] for p in pairs}:
                fs.write(f">{k}\n{sseqs[k]}\n")
        out = subprocess.run(minimap2.split() + ["-c", "--secondary=yes", "-N", "200", "-t", str(threads), sf, qf],
                             check=True, capture_output=True, text=True).stdout
        for line in out.splitlines():
            f = line.split("\t")
            key = (f[0], f[5])
            if key not in pairs or f[4] != "+":
                continue
            qlen, tlen, nm, al = int(f[1]), int(f[6]), int(f[9]), int(f[10])
            cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
            aligned = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg) if op in "M=X")   # target bases
            if key not in sim or nm > sim[key][3]:                                                  # aligned, no gaps
                sim[key] = (aligned / tlen, qlen / tlen, 100.0 * nm / al, nm)
        sim = {k: v[:3] for k, v in sim.items()}
        rest = {p for p in pairs if p not in sim or min(len(qseqs[p[0]]), len(sseqs[p[1]])) < short_cds}
        if rest:
            with open(qf, "w") as fq:
                for k in {p[0] for p in rest}:
                    fq.write(f">{k}\n{qseqs[k]}\n")
            with open(sf, "w") as fs:
                for k in {p[1] for p in rest}:
                    fs.write(f">{k}\n{sseqs[k]}\n")
            subprocess.run(["makeblastdb", "-in", sf, "-dbtype", "nucl", "-out", os.path.join(tmp, "db")], check=True,
                           stdout=subprocess.DEVNULL)
            out = subprocess.run(["blastn", "-query", qf, "-db", os.path.join(tmp, "db"), "-evalue", "1e-10",
                                  "-dust", "no", "-max_target_seqs", "500", "-num_threads", str(threads), "-outfmt",
                                  "6 qseqid sseqid pident length qlen slen qstart qend sstart send bitscore"],
                                 check=True, capture_output=True, text=True).stdout
            best = {}
            for line in out.splitlines():
                f = line.split("\t")
                key = (f[0], f[1])
                if key not in rest or int(f[8]) > int(f[9]):
                    continue
                if key not in best or float(f[10]) > best[key][0]:
                    best[key] = (float(f[10]), (int(f[9]) - int(f[8]) + 1) / int(f[5]), int(f[4]) / int(f[5]),
                                 float(f[2]))
            for k, v in best.items():
                sim[k] = v[1:]
    return sim


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baskets", required=True, help="baskets_v2.tsv (winners, n_members, members)")
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--transfer-dir", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--max-len", type=float, default=1.2)
    ap.add_argument("--min-id", type=float, default=95.0, help="keep_similar: minimum identity of the own gene")
    ap.add_argument("--keep-min-cov", type=float, default=0.8,
                    help="keep_similar: the own gene must cover >= this share of the target CDS")
    ap.add_argument("--min-cov", "--min-cov-target", dest="min_cov", type=float, default=0.7,
                    help="a valid candidate must cover >= this share of the target CDS")
    ap.add_argument("--min-cov-rescue", type=float, default=0.5, help="... or >= this for a rescue")
    ap.add_argument("--min-tcov", type=float, default=0.7,
                    help="the target's transfer without a valid ORF is kept (target_invalid) when its coverage is >= this")
    ap.add_argument("--min-tcov-rescue", type=float, default=0.5,
                    help="another member's transfer is a candidate when its Liftoff coverage is >= this")
    ap.add_argument("--contig-end", type=int, default=1000, help="a model this close to a contig end counts as valid")
    ap.add_argument("--short-cds", type=int, default=300, help="CDS shorter than this are compared with blastn")
    ap.add_argument("--minimap2", default="minimap2 -x asm10 -k 13 -w 5 -s 20 -m 20", help="minimap2 command (settings for short CDS)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--arrays", default="", help="baskets_v2.py's tandem_arrays.tsv: baskets holding a tandem array")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]

    baskets = []
    for r in csv.DictReader(open(a.baskets), delimiter="\t"):
        baskets.append(([w for w in r["winners"].split(";") if w], r["members"].split(",")))   # no winner: ""
    allm = {m for _, ms in baskets for m in ms}
    arrays = {}
    if a.arrays and os.path.exists(a.arrays):
        arrays = {r["winners"]: set(r["copies"].split(",")) for r in csv.DictReader(open(a.arrays), delimiter="\t")}
    winners = {w for ws, _ in baskets for w in ws}

    native, owner, seqs = {}, {}, {}
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))
        for k, v in read_models(path, allm).items():
            native[k] = v
            owner[k] = g
    print(f"{len(native)} native models", file=sys.stderr)
    tr = {}
    for path in sorted(glob.glob(os.path.join(a.transfer_dir, "*_onto_*"))):
        src, tgt = os.path.basename(path).split("_liftoff")[0].split("_onto_")
        if tgt not in genomes:
            continue
        want = {m for m in allm if owner.get(m) == src}
        for k, v in read_models(path, want, primary_only=True).items():
            tr[(k, tgt)] = v
    print(f"{len(tr)} transfers", file=sys.stderr)

    # sequences and ORFs, one genome at a time (st: valid / gap / end / invalid)
    orf, st = {}, {}
    for g in genomes:
        genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
        clen = {c: len(s) for c, s in genome.items()}
        for k, m in native.items():
            if owner[k] == g:
                s = spliced(genome, m)
                seqs[("own", k)] = s
                st[("own", k)] = orf_status(s, m, clen, a.contig_end)
                orf[("own", k)] = st[("own", k)] == "valid"
        for (k, t), m in tr.items():
            if t == g:
                st[("tr", k, t)] = orf_status(spliced(genome, m), m, clen, a.contig_end)
                orf[("tr", k, t)] = st[("tr", k, t)] == "valid"
        del genome
        print(f"  ORFs of {g}: done", file=sys.stderr)

    # every genome's own member CDS against the target CDS of its basket (for keep_similar)
    pairs = set()
    for ws, ms in baskets:
        for w in ws:
            for m in ms:
                if m != w and owner.get(m) and owner.get(w):
                    pairs.add((m, w))
    wseqs = {w: seqs.get(("own", w)) for w in winners}
    sim = cds_sim({m: seqs.get(("own", m)) for m, _ in pairs}, wseqs, pairs, a.threads, a.minimap2, a.short_cds)
    print(f"{len(sim)} of {len(pairs)} member-target pairs aligned", file=sys.stderr)

    claims = Counter()   # (own gene) -> number of panel genes whose target model, in its genome, lies over it
    for ws, ms in baskets:
        for w in ws:
            for m in ms:
                t = owner.get(m)
                if not t or m == w or owner.get(w) == t:
                    continue
                tm = tr.get((w, t))
                if tm and ov(span(native[m]), span(tm)) > 0:
                    claims[m] += 1

    # candidates of each (panel gene, genome): the target's transfer, the genome's own members, the other members'
    # transfers (Liftoff coverage >= --min-tcov-rescue), at the target's place. Absent cases keep why for absent.tsv.
    cands, why_none = defaultdict(list), {}
    for ws, ms in baskets:
        for w in ws:
            for t in genomes:
                if owner[w] == t:
                    continue
                tm = tr.get((w, t))
                anchors = [span(native[m]) for m in ms if owner.get(m) == t and m != w]
                n_tr = n_cov = n_place = 0
                if tm is not None:
                    n_tr += 1
                    n_cov += 1
                    n_place += 1
                    cands[(w, t)].append(("tr", w))
                for m in ms:
                    if m == w:
                        continue
                    if owner.get(m) == t:
                        if claims[m] > 1:
                            continue
                        mm, key = native[m], ("own", m)
                    else:
                        mm, key = tr.get((m, t)), ("tr", m)
                        if mm is None:
                            continue
                        n_tr += 1
                        if mm[3] < a.min_tcov_rescue:
                            continue
                        n_cov += 1
                    sm = span(mm)
                    if tm is not None:
                        if ov(sm, span(tm)) == 0:
                            continue
                    elif key[0] == "tr" and anchors and not any(ov(sm, x) > 0 for x in anchors):
                        continue
                    n_place += 1
                    cands[(w, t)].append(key)
                why_none[(w, t)] = (n_tr, n_cov, n_place)
    # sequences of the candidate transfers, then their similarity to the target CDS
    want = defaultdict(set)
    for (w, t), ks in cands.items():
        for kind, m in ks:
            if kind == "tr":
                want[t].add(m)
    cseq = {}
    for g in genomes:
        if not want[g]:
            continue
        genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
        for m in want[g]:
            cseq[f"{m}@{g}"] = spliced(genome, tr[(m, g)])
        del genome
    cpairs, qs = set(), {}
    for (w, t), ks in cands.items():
        for kind, m in ks:
            q = m if kind == "own" else f"{m}@{t}"
            qs[q] = seqs.get(("own", m)) if kind == "own" else cseq.get(q)
            cpairs.add((q, w))
    csim = cds_sim(qs, wseqs, cpairs, a.threads, a.minimap2, a.short_cds)
    print(f"{sum(len(v) for v in cands.values())} candidate models for {len(cands)} (panel gene, genome) pairs; "
          f"{len(csim)} aligned", file=sys.stderr)

    def scored(w, t):
        """the candidates of (w, t): (key, kind, member, ORF status, coverage, length ratio, identity)."""
        out = []
        for kind, m in cands.get((w, t), []):
            q = m if kind == "own" else f"{m}@{t}"
            s_ = st.get(("own", m)) if kind == "own" else st.get(("tr", m, t))
            sc = csim.get((q, w))
            if sc:
                out.append((kind, m, s_, sc[0], sc[1], sc[2]))
        return out

    def best_of(w, cs, min_cov):
        """the most similar valid candidate covering >= min_cov of the target CDS and <= --max-len of its length;
        on a tie, the target's own transfer."""
        best, key = None, None
        for c in cs:
            kind, m, s_, cov, ratio, ident = c
            if s_ not in VALID_LIKE or cov < min_cov or ratio > a.max_len:
                continue
            k = (round(cov, 3), round(ident, 2), -round(abs(1 - ratio), 3), kind == "tr" and m == w)
            if key is None or k > key:
                best, key = c, k
        return best

    # tandem arrays: besides the normal model of each winner, a genome keeps its own confirmed copies (baskets_v2.py's
    # copies), and gets the transfers of other genomes' copies that have a valid ORF and land where it has no gene of
    # the basket (Liftoff coverage >= --min-tcov; of transfers overlapping each other, the highest coverage, then the
    # longest). The build adds them where they do not overlap the models already written for the basket.
    array_rows = []
    for ws, ms in baskets:
        key = ";".join(ws)
        if key not in arrays:
            continue
        cp = arrays[key]
        for t in genomes:
            own = [m for m in ms if owner.get(m) == t]
            array_rows += [[key, t, m, "own", st.get(("own", m), "")] for m in own if m in cp]
            fills = []
            for m in cp:
                mm = tr.get((m, t))
                if owner.get(m) == t or mm is None or mm[3] < a.min_tcov or st.get(("tr", m, t)) != "valid":
                    continue
                if any(ov(span(mm), span(native[k])) > 0 for k in own):
                    continue
                fills.append(m)
            fills.sort(key=lambda m: (tr[(m, t)][3], sum(b - a_ + 1 for a_, b in tr[(m, t)][2])), reverse=True)
            chosen = []
            for m in fills:
                if not any(ov(span(tr[(m, t)]), span(tr[(c, t)])) > 0 for c in chosen):
                    chosen.append(m)
            array_rows += [[key, t, m, "transfer", st.get(("tr", m, t), "")] for m in chosen]
    with open(os.path.join(a.out_dir, "array_models.tsv"), "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["winners", "genome", "member", "model", "orf_status"])
        wr.writerows(array_rows)
    if arrays:
        print(f"{len(arrays)} tandem-array baskets: {sum(1 for r in array_rows if r[3] == 'own')} own copies, "
              f"{sum(1 for r in array_rows if r[3] == 'transfer')} valid transfers of copies, to add where they fit",
              file=sys.stderr)

    rows, summ, absent_rows = [], defaultdict(Counter), []
    for ws, ms in baskets:
        for w in ws:
            for t in genomes:
                if owner[w] == t:
                    tm, tv = native[w], orf[("own", w)]
                else:
                    tm = tr.get((w, t))
                    tv = orf.get(("tr", w, t)) if tm else None
                own = [m for m in ms if owner.get(m) == t and m != w and
                       (tm is None or ov(span(native[m]), span(tm)) > 0)]
                if owner[w] == t:
                    own = []
                alt, why, chosen = "", "", ""
                if owner[w] == t:
                    d_panel = "target (own model)"   # a gene is never replaced in its own genome, even with an invalid ORF
                    chosen = st.get(("own", w), "")
                else:
                    cs = scored(w, t)
                    b = best_of(w, cs, a.min_cov)
                    if b and b[0] == "tr" and b[1] == w:
                        d_panel, chosen = ("replace" if own else "fill"), b[2]
                    elif b:
                        d_panel, alt, chosen = "alternative", b[1], b[2]
                    elif tm is not None and tm[3] >= a.min_tcov and st.get(("tr", w, t)) != "valid":
                        d_panel, chosen = "target_invalid", st.get(("tr", w, t), "")
                    else:
                        b = best_of(w, cs, a.min_cov_rescue)
                        if b:
                            d_panel, alt, chosen = "rescued", b[1], b[2]
                        elif own:
                            d_panel = "keep_divergent"
                            chosen = st.get(("own", max(own, key=lambda k: span(native[k])[3] - span(native[k])[2])), "")
                        else:
                            d_panel = "absent"
                            n_tr, n_cov, n_place = why_none.get((w, t), (0, 0, 0))
                            okc = [c for c in cs if c[2] in VALID_LIKE]
                            if n_tr == 0:
                                reason = "no transfer of any member into this genome"
                            elif n_cov == 0:
                                reason = f"transfers with Liftoff coverage < {a.min_tcov_rescue}"
                            elif n_place == 0:
                                reason = "transfers not at the panel gene's place"
                            elif not cs:
                                reason = "models at the place do not align with the panel gene's CDS"
                            elif not okc:
                                reason = "models at the place without a valid ORF"
                            else:
                                reason = f"valid models cover < {a.min_cov_rescue} of the panel gene's CDS (or too long)"
                            top = max(cs, key=lambda c: (c[2] in VALID_LIKE, c[3]), default=None)
                            absent_rows.append([w, t, reason, f"{tm[3]:.3f}" if tm else "",
                                                st.get(("tr", w, t), "") if tm else "",
                                                (top[1] if top[0] == "own" else f"{top[1]} (transfer)") if top else "",
                                                f"{top[3]:.3f}" if top else "", top[2] if top else ""])
                # mode "keep_similar"
                d_keep = d_panel
                if len(own) == 1 and owner[w] != t and d_panel not in ("keep_divergent", "absent") and claims[own[0]] > 1:
                    why = "own gene lies under 2+ panel genes (a fusion): not kept"
                elif len(own) == 1 and owner[w] != t and d_panel not in ("keep_divergent", "absent"):
                    m = own[0]
                    sc = sim.get((m, w))
                    if sc:
                        why = f"cov {sc[0]:.2f} len {sc[1]:.2f} id {sc[2]:.1f}"
                    if orf.get(("own", m)) and sc and sc[0] >= a.keep_min_cov and sc[1] <= a.max_len and sc[2] >= a.min_id:
                        d_keep = "keep_original"
                elif len(own) > 1:
                    why = f"{len(own)} own genes (pieces or copies)"
                summ[t][("prefer_panel_model", d_panel)] += 1
                summ[t][("keep_similar", d_keep)] += 1
                rows.append([w, t, ",".join(own), d_panel, d_keep, why, alt, chosen])
    with open(os.path.join(a.out_dir, "absent.tsv"), "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["panel_gene", "genome", "reason", "panel_gene_transfer_coverage", "panel_gene_transfer_orf",
                     "best_candidate", "best_candidate_cov_of_panel_cds", "best_candidate_orf"])
        wr.writerows(absent_rows)
    with open(os.path.join(a.out_dir, "own_invalid_orf.tsv"), "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["gene", "genome", "is_panel_gene", "start_ok", "stop_ok", "frame_ok", "internal_stops"])
        for k in sorted(native):
            sq = seqs.get(("own", k))
            if sq and not orf.get(("own", k)):
                st_, sp_, fr_, inn_, _ = check_orf.check(sq)
                wr.writerow([k, owner[k], int(k in winners), st_, sp_, fr_, inn_])
    with open(os.path.join(a.out_dir, "decisions.tsv"), "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["panel_gene", "genome", "own_genes", "decision_prefer_panel_model", "decision_keep_similar", "comparison",
                     "alternative_member", "model_orf_prefer_panel_model"])
        wr.writerows(rows)
    keys = sorted({k for c in summ.values() for k in c})
    with open(os.path.join(a.out_dir, "summary.tsv"), "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["genome"] + [f"{m}:{d}" for m, d in keys])
        for t in genomes:
            wr.writerow([t] + [summ[t][k] for k in keys])
    tot = Counter()
    for c in summ.values():
        tot.update(c)
    for mode in ("keep_similar", "prefer_panel_model"):
        print(f"\n== mode {mode}")
        for (m, d), v in sorted(tot.items()):
            if m == mode:
                print(f"  {v:7}  {d}")


if __name__ == "__main__":
    main()
