#!/usr/bin/env python3
"""
Classify blastn hits of the reference organellar sequence(s) against each genome.

Inputs:
  --blast-dir   one <label>.tsv per genome, blastn -outfmt
                '6 qseqid sseqid pident length mismatch gapopen qstart qend
                    sstart send evalue bitscore qlen slen'
  --organelles  the organelles.fna those searches used (for the organelle lengths)

Output: one row per (genome, organelle, contig) that has a real hit, so a genome
whose assembly has the organelle split — or duplicated — over several contigs
shows up as several rows.

  contig_class  full              contig ~= one organelle copy (>= --min-cov of
                                  the organelle, contig <= 1.5x its length)
                full_in_scaffold  whole organelle, but embedded in a large contig
                fragment          partial hit (NUMT / NUPT, or a broken piece)

  genome_call   own_contig            exactly one 'full' contig
                MULTIPLE_COPIES       >1 'full' contig  <- usually an assembly
                                      artefact; check before keeping any of them
                split_across_contigs  no single full contig, but the pieces
                                      together cover >= --min-cov
                in_scaffold           a 'full_in_scaffold' and no 'full'
                fragments_only        only fragments, together < --min-cov
                absent                no hit at all

Also prints a per-organelle summary.
"""
import argparse
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.intervals import merged_len


def organelle_lengths(path):
    lengths, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast-dir", required=True)
    ap.add_argument("--organelles", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cov", type=float, default=0.90,
                    help="organelle-coverage fraction for a 'full' contig (default 0.90)")
    ap.add_argument("--min-partial", type=float, default=0.05,
                    help="coverage below which a contig hit is ignored (default 0.05)")
    a = ap.parse_args()

    olen = organelle_lengths(a.organelles)
    if not olen:
        raise SystemExit(f"no sequences in {a.organelles}")

    # (genome, organelle, contig) -> query intervals on the organelle, id-weighted
    #                                sums, contig length
    qiv = defaultdict(list)
    id_wsum = defaultdict(float)
    id_wlen = defaultdict(float)
    contig_len = {}
    genomes = set()

    for tsv in sorted(glob.glob(os.path.join(a.blast_dir, "*.tsv"))):
        genome = os.path.basename(tsv)[:-4]
        genomes.add(genome)
        with open(tsv) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 14:
                    continue
                organelle, contig = f[0], f[1]
                ident = float(f[2])
                qs, qe = sorted((int(f[6]), int(f[7])))
                k = (genome, organelle, contig)
                qiv[k].append((qs, qe))
                id_wsum[k] += ident * (qe - qs)
                id_wlen[k] += (qe - qs)
                contig_len[k] = int(f[13])

    # per-contig stats + class
    contig_stat = {}
    for k, ivs in qiv.items():
        _, organelle, _ = k
        ol = olen.get(organelle) or 1
        cov = merged_len(ivs) / ol
        pid = id_wsum[k] / id_wlen[k] if id_wlen[k] else 0.0
        clen = contig_len[k]
        if cov >= a.min_cov and clen <= 1.5 * ol:
            cls = "full"
        elif cov >= a.min_cov:
            cls = "full_in_scaffold"
        elif cov >= a.min_partial:
            cls = "fragment"
        else:
            cls = "trace"
        contig_stat[k] = dict(cov=cov, pid=pid, clen=clen, ol=ol, cls=cls)

    # per (genome, organelle) call
    per_go = defaultdict(list)
    for k, v in contig_stat.items():
        per_go[(k[0], k[1])].append((k[2], v))

    genome_call = {}
    for (g, o), lst in per_go.items():
        ol = olen.get(o) or 1
        fulls = [c for c, v in lst if v["cls"] == "full"]
        scaf = [c for c, v in lst if v["cls"] == "full_in_scaffold"]
        union = merged_len([iv for k, ivs in qiv.items()
                            if k[0] == g and k[1] == o for iv in ivs]) / ol
        if len(fulls) > 1:
            call = "MULTIPLE_COPIES"
        elif len(fulls) == 1:
            call = "own_contig"
        elif scaf:
            call = "in_scaffold"
        elif union >= a.min_cov:
            call = "split_across_contigs"
        elif any(v["cls"] == "fragment" for _, v in lst):
            call = "fragments_only"
        else:
            call = "absent"
        genome_call[(g, o)] = call

    rows = []
    for k, v in contig_stat.items():
        if v["cls"] == "trace":
            continue
        g, o, c = k
        rows.append([g, o, c, round(v["pid"], 1), round(v["cov"], 3),
                     v["clen"], v["ol"], v["cls"], genome_call[(g, o)]])
    for g in genomes:
        for o in olen:
            if (g, o) not in genome_call:
                rows.append([g, o, "-", 0.0, 0.0, 0, olen[o], "-", "absent"])
                genome_call[(g, o)] = "absent"

    rows.sort(key=lambda r: (r[1], r[0], -r[4]))
    with open(a.out, "w") as out:
        out.write("genome\torganelle\tcontig\tpident\tcov_of_organelle\t"
                  "contig_len\torganelle_len\tcontig_class\tgenome_call\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {a.out}\n")

    for o in olen:
        gc = {g: genome_call[(g, o)] for g in genomes}
        def names(*calls):
            return sorted(g for g, c in gc.items() if c in calls) or ["-"]
        multi = sorted(
            (g, sum(1 for _, v in per_go[(g, o)] if v["cls"] == "full"))
            for g, c in gc.items() if c == "MULTIPLE_COPIES")
        print(f"[{o}]  (len {olen[o]})")
        print(f"   own contig (1 copy)  : {', '.join(names('own_contig'))}")
        print(f"   MULTIPLE copies      : "
              f"{', '.join(f'{g} ({n} contigs)' for g, n in multi) or '-'}   <- check")
        print(f"   inside a scaffold    : {', '.join(names('in_scaffold'))}")
        print(f"   split across contigs : {', '.join(names('split_across_contigs'))}")
        print(f"   only fragments       : {', '.join(names('fragments_only'))}")
        print(f"   no hit               : {', '.join(names('absent'))}\n")


if __name__ == "__main__":
    main()
