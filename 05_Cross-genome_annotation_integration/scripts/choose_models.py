#!/usr/bin/env python3
"""
Test of steps 4-5 of REVISION_PLAN.md, on the baskets of build_baskets.py. Counts only; no GFF is written.

Step 4, one model per basket, applied to every genome:
- basket with a reference gene R: genome T takes R's transfer when it has a valid ORF and coverage >= --min-cov;
  otherwise the best member whose transfer onto T has a valid ORF and overlaps R's transfer by >= --maps-frac of each
  (it "maps to R"); otherwise R's transfer, marked ref_invalid_orf. A basket is "ref model doubtful" when R's transfer
  has no valid ORF in most genomes where it lands, while another member has a valid ORF in most genomes.
- basket without a reference gene: the member whose transfer has a valid ORF (coverage >= --min-cov) in the most
  genomes wins (then mean identity, then CDS span); each genome takes its transfer, or the best valid alternative.
- a genome gets the gene only if the chosen transfer covers >= --min-cov; in its own genome a member is its native
  model (no ORF check here: counted as valid when its transfers are valid in most genomes).
Step 5: pairs of baskets whose chosen models overlap (CDS spans, >= 20 bp) in each genome, and in how many genomes.

Usage: python3 choose_models.py --baskets baskets.tsv --native-dir gff_originals --transfer-dir polished_nonredundant \
           --genomes genome_order.txt --out-dir test
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_baskets import read_gff   # noqa: E402

VO = re.compile(r"valid_ORFs?=(\d+)")
COV = re.compile(r"(?:^|;)coverage=([\d.]+)")
IDN = re.compile(r"(?:^|;)sequence_ID=([\d.]+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baskets", required=True)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--transfer-dir", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--reference", default="")
    ap.add_argument("--min-cov", type=float, default=0.8)
    ap.add_argument("--maps-frac", type=float, default=0.8)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    ref = a.reference or genomes[0]

    baskets, piece_of = {}, {}
    for r in csv.DictReader(open(a.baskets), delimiter="\t"):
        baskets[r["basket"]] = r["natives"].split(",")
        for x in r["partial_members"].split(","):
            if x:
                piece_of.setdefault(x, r["basket"])
    # a basket made only of genes that are pieces of another basket is not a gene of its own: it is left out
    pieces_only = {b for b, mem in baskets.items() if all(m in piece_of and piece_of[m] != b for m in mem)}
    print(f"baskets made only of pieces of another basket (left out): {len(pieces_only)} of {len(baskets)}")
    baskets = {b: m for b, m in baskets.items() if b not in pieces_only}
    owner, native = {}, {}
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))
        for gid, (c, st, s, e, _) in read_gff(path).items():
            owner[gid] = g
            native[gid] = (c, st, s, e)
    # transfers: (source gene, target) -> (contig, strand, s, e, valid, cov, id)
    tr = {}
    for path in sorted(glob.glob(os.path.join(a.transfer_dir, "*_onto_*"))):
        src, tgt = os.path.basename(path).split("_liftoff")[0].split("_onto_")
        if tgt not in genomes:
            continue
        for gid, (c, st, s, e, attr) in read_gff(path, primary_only=True).items():
            m, cv, idn = VO.search(attr), COV.search(attr), IDN.search(attr)
            tr[(gid, tgt)] = (c, st, s, e, (m.group(1) != "0") if m else False,
                              float(cv.group(1)) if cv else 0.0, float(idn.group(1)) if idn else 0.0)
    print(f"{len(tr)} transfers", file=sys.stderr)

    def model(s, t):
        """the model member s gives genome t: (contig, strand, s, e, valid, cov, id) or None."""
        if owner[s] == t:
            c, st, x, y = native[s]
            v = [tr[(s, h)][4] for h in genomes if (s, h) in tr]
            return (c, st, x, y, (sum(v) * 2 > len(v)) if v else True, 1.0, 1.0)
        return tr.get((s, t))

    def ov(m1, m2):
        if m1[0] != m2[0] or m1[1] != m2[1]:
            return 0, 0, 0
        o = min(m1[3], m2[3]) - max(m1[2], m2[2]) + 1
        return o, m1[3] - m1[2] + 1, m2[3] - m2[2] + 1

    stats = Counter()
    chosen_models = defaultdict(list)   # genome -> (contig, strand, s, e, basket)
    rows = []
    for b, mem in baskets.items():
        refs = [m for m in mem if owner.get(m) == ref]
        has_ref = bool(refs)
        if has_ref:
            R = refs[0]
            panel = R
        else:
            def score(s):
                ms = [model(s, t) for t in genomes]
                ok = [m for m in ms if m and m[4] and m[5] >= a.min_cov]
                ids = [m[6] for m in ms if m]
                return (len(ok), sum(ids) / len(ids) if ids else 0, native[s][3] - native[s][2])
            panel = max(mem, key=score)
        states = Counter()
        r_invalid = r_present = 0
        for t in genomes:
            m = model(panel, t)
            if has_ref and m and m[5] >= a.min_cov:
                r_present += 1
                r_invalid += not m[4]
            if m and m[5] >= a.min_cov and m[4]:
                use, why = m, "panel model"
            else:
                alts = []
                for s in mem:
                    if s == panel:
                        continue
                    ms = model(s, t)
                    if not ms or not ms[4] or ms[5] < a.min_cov:
                        continue
                    if has_ref and m and m[5] >= a.min_cov:
                        o, l1, l2 = ov(ms, m)
                        if o < a.maps_frac * l1 or o < a.maps_frac * l2:
                            continue
                    alts.append((ms[5] * ms[6], ms))
                if alts:
                    use, why = max(alts, key=lambda x: x[0])[1], "replacement with a valid ORF"
                elif m and m[5] >= a.min_cov:
                    use, why = m, ("ref_invalid_orf" if has_ref else "panel model, invalid ORF")
                else:
                    use, why = None, "absent"
            states[why] += 1
            if use:
                chosen_models[t].append((use[0], use[1], use[2], use[3], b))
        doubtful = False
        if has_ref and r_present and r_invalid * 2 > r_present:
            for s in mem:
                if s == R:
                    continue
                v = [model(s, t) for t in genomes]
                v = [x for x in v if x and x[5] >= a.min_cov]
                if v and sum(x[4] for x in v) * 2 > len(v):
                    doubtful = True
                    break
        kind = "with reference" if has_ref else ("no reference, 1 native" if len(mem) == 1 else "no reference, 2+ natives")
        stats[(kind, "baskets")] += 1
        for k, v in states.items():
            stats[(kind, k)] += v
        stats[(kind, "ref model doubtful")] += doubtful
        present = sum(v for k, v in states.items() if k != "absent")
        stats[(kind, f"genomes with the gene: {'15' if present == len(genomes) else '2-14' if present > 1 else present}")] += 1
        rows.append([b, kind, panel, present, states["panel model"], states["replacement with a valid ORF"],
                     states["ref_invalid_orf"] + states["panel model, invalid ORF"], states["absent"], int(doubtful)])
    with open(os.path.join(a.out_dir, "choices.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["basket", "kind", "panel_model", "genomes_with_gene", "panel_model_used", "replacement",
                    "invalid_kept", "absent", "ref_model_doubtful"])
        w.writerows(rows)
    print("\n== step 4: one model per basket")
    for kind in ("with reference", "no reference, 2+ natives", "no reference, 1 native"):
        print(f"-- {kind}")
        for k in sorted(k2 for (k1, k2) in stats if k1 == kind):
            print(f"   {stats[(kind, k)]:8}  {k}")

    # step 5: overlaps between baskets
    pairs = defaultdict(set)
    for t, L in chosen_models.items():
        L.sort(key=lambda x: (x[0], x[2]))
        act = []
        for c, st, s, e, b in L:
            act = [x for x in act if x[0] == c and x[3] >= s]
            for c2, st2, s2, e2, b2 in act:
                if b2 != b and min(e, e2) - max(s, s2) + 1 >= 20:
                    inside = (s >= s2 and e <= e2) or (s2 >= s and e2 <= e)
                    pairs[(min(b, b2), max(b, b2), "same" if st == st2 else "opposite",
                           "contained" if inside else "partial")].add(t)
            act.append((c, st, s, e, b))
    print("\n== step 5: baskets whose chosen models overlap")
    by = Counter()
    for (b1, b2, strand, kind), gs in pairs.items():
        n = len(gs)
        by[(strand, kind, "1-2 genomes" if n <= 2 else "3-7" if n <= 7 else "8-15")] += 1
    for k in sorted(by):
        print(f"   {by[k]:6}  {k[0]} strand, {k[1]}, overlapping in {k[2]}")
    print(f"   total basket pairs: {len(pairs)}")
    with open(os.path.join(a.out_dir, "basket_overlaps.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["basket1", "basket2", "strand", "type", "n_genomes", "genomes"])
        for (b1, b2, strand, kind), gs in sorted(pairs.items(), key=lambda x: -len(x[1])):
            w.writerow([b1, b2, strand, kind, len(gs), ",".join(sorted(gs))])


if __name__ == "__main__":
    main()
