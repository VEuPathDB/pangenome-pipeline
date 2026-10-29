#!/usr/bin/env python3
"""Check of the tandem-array detection of baskets_v2.py: the same neighbouring pairs (two genes of one genome in one
basket, same strand, not overlapping, at most --array-k genes apart, protein lengths within --array-cov), their blastp
identity and coverage as baskets_v2.py computes them, and, for the --recheck pairs closest to the thresholds, an
independent blastp of the two proteins alone (-subject), to confirm the scores are read correctly.

Writes <out>.tsv (one row per pair) and prints a summary."""
import argparse
import glob
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_baskets import read_gff   # noqa: E402
from decide_models import read_models, spliced   # noqa: E402
from baskets_v2 import protein, protein_pairs   # noqa: E402
import check_orf   # noqa: E402


def blastp_alone(px, py):
    """best HSP of px against py, blastp -subject: (identity %, share of x aligned, share of y aligned)."""
    with tempfile.TemporaryDirectory() as tmp:
        fx, fy = os.path.join(tmp, "x.faa"), os.path.join(tmp, "y.faa")
        open(fx, "w").write(f">x\n{px}\n")
        open(fy, "w").write(f">y\n{py}\n")
        res = subprocess.run(["blastp", "-query", fx, "-subject", fy, "-evalue", "1e-5", "-outfmt",
                              "6 pident qstart qend sstart send qlen slen bitscore"],
                             check=True, capture_output=True, text=True).stdout
    best = None
    for line in res.splitlines():
        f = line.split("\t")
        if best is None or float(f[7]) > best[0]:
            best = (float(f[7]), float(f[0]), (int(f[2]) - int(f[1]) + 1) / int(f[5]),
                    (int(f[4]) - int(f[3]) + 1) / int(f[6]))
    return best[1:] if best else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baskets", required=True, help="baskets_v2.tsv")
    ap.add_argument("--native-dir", required=True, help="the native GFFs given to baskets_v2.py")
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", default=".fna")
    ap.add_argument("--array-k", type=int, default=3)
    ap.add_argument("--array-id", type=float, default=70.0)
    ap.add_argument("--array-cov", type=float, default=0.7)
    ap.add_argument("--recheck", type=int, default=30, help="pairs rechecked with a blastp of the two proteins alone")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out", default="tandem_pairs_check")
    a = ap.parse_args()

    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    ref = genomes[0]
    native, owner, nmodel = {}, {}, {}
    for g in genomes:
        path = next(iter(glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff"))))
        genes = read_gff(path)
        for gid, (c, st, s, e, _) in genes.items():
            native[gid] = (c, st, s, e)
            owner[gid] = g
        nmodel.update(read_models(path, set(genes)))
    order = defaultdict(list)
    for gid, (c, st, s, e) in native.items():
        order[(owner[gid], c)].append((s, gid))
    pos = {}
    for k in order:
        order[k].sort()
        for i, (_, gid) in enumerate(order[k]):
            pos[gid] = (k, i)

    baskets = []
    for line in list(open(a.baskets))[1:]:
        w, _, ms = line.rstrip("\n").split("\t")
        if w:
            baskets.append((w, [m for m in ms.split(",") if m in native]))
    need = {m for _, ms in baskets for m in ms if owner[m] != ref}
    prot = {}
    for g in genomes:
        if g == ref:
            continue
        genome = check_orf.read_fasta(os.path.join(a.genome_dir, g + a.genome_suffix))
        for gid in need:
            if owner[gid] == g and gid in nmodel:
                prot[gid] = protein(spliced(genome, nmodel[gid]))
        del genome
        print(f"  proteins of {g}: done", file=sys.stderr, flush=True)

    pairs = []
    for key, ms in baskets:
        by_genome = defaultdict(list)
        for m in ms:
            by_genome[owner[m]].append(m)
        for g, gm in by_genome.items():
            for i, x in enumerate(gm):
                for y in gm[i + 1:]:
                    (kx, ix), (ky, iy) = pos[x], pos[y]
                    nx, ny = native[x], native[y]
                    if kx != ky or nx[1] != ny[1] or min(nx[3], ny[3]) >= max(nx[2], ny[2]):
                        continue
                    lx, ly = len(prot.get(x, "")), len(prot.get(y, ""))
                    if abs(ix - iy) <= a.array_k and lx and ly and min(lx, ly) >= a.array_cov * max(lx, ly):
                        pairs.append((key, g, x, y))
    by_genome = defaultdict(set)
    for key, g, x, y in pairs:
        by_genome[g].add((x, y))
    print(f"{len(pairs)} neighbouring pairs, {len({m for v in by_genome.values() for p in v for m in p})} proteins",
          file=sys.stderr, flush=True)
    sim = protein_pairs(prot, by_genome, a.threads)

    rows = []
    for key, g, x, y in pairs:
        sc = sim.get((x, y)) or sim.get((y, x))
        rows.append((key, g, x, y, len(prot[x]), len(prot[y])) + (sc if sc else (0.0, 0.0, 0.0)))
    rows.sort(key=lambda r: (min(r[6] / a.array_id, r[7] / a.array_cov, r[8] / a.array_cov), r[6]), reverse=True)

    passed = Counter()
    with open(a.out + ".tsv", "w") as fh:
        fh.write("winners\tgenome\tgene_x\tgene_y\taa_x\taa_y\tidentity\tcov_x\tcov_y\tpasses\trecheck_identity\t"
                 "recheck_cov_x\trecheck_cov_y\n")
        for n, r in enumerate(rows):
            ok = r[6] >= a.array_id and r[7] >= a.array_cov and r[8] >= a.array_cov
            passed["passes" if ok else "no blastp hit" if r[6] == 0 else
                   "identity below" if r[6] < a.array_id else "coverage below"] += 1
            rc = ""
            if n < a.recheck:
                b = blastp_alone(prot[r[2]], prot[r[3]])
                rc = "\t".join(f"{v:.3f}" for v in b) if b else "no hit\t\t"
            fh.write("\t".join(map(str, r[:6])) + f"\t{r[6]:.1f}\t{r[7]:.3f}\t{r[8]:.3f}\t{'yes' if ok else 'no'}\t"
                     + (rc or "\t\t") + "\n")

    print(f"{len(rows)} pairs: " + ", ".join(f"{k} {v}" for k, v in passed.most_common()))
    ids = sorted(r[6] for r in rows if r[6] > 0)
    for t in (90, 80, 70, 60, 50, 40):
        print(f"  identity >= {t}%: {sum(1 for v in ids if v >= t)} pairs "
              f"({sum(1 for r in rows if r[6] >= t and min(r[7], r[8]) >= a.array_cov)} also covering >= "
              f"{a.array_cov:.0%} of both)")
    print(f"closest pairs (and their recheck) in {a.out}.tsv")


if __name__ == "__main__":
    main()
