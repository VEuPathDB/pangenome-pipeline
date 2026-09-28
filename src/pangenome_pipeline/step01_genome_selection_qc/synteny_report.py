#!/usr/bin/env python3
"""
Summarise minimap2 PAF whole-genome alignments vs a reference, one row per query
genome, with orientation-corrected rearrangement detection.

    synteny_report.py <synteny_dir> [--min-aln 2000] [--ref-mb 65.63]
                      [--md out.md] [--tsv out.tsv]

<synteny_dir> holds  <label>_vs_<REF>.paf  files (as written by synteny_vs_ref.slurm.sh).

Why orientation correction: plot_dotplot.py colours a whole query scaffold red
whenever that scaffold was assembled opposite to the reference. For fragmented
(scaffold-level) assemblies ~half of all scaffolds flip by chance, so raw
"% reverse-strand" is meaningless. Here each query contig is first assigned its
own dominant orientation (by aligned bp to its main reference chromosome); only
alignment blocks that run *against their own contig's* orientation, on the same
reference chromosome, are counted as candidate internal inversions.
"""

import argparse
import collections
import glob
import os
import re
import sys

from pangenome_pipeline.common.intervals import merged_len
from pangenome_pipeline.common.io import read_paf

TRANSLOC_MIN = 20_000  # bp to a 2nd ref chrom before a contig counts as split/chimeric
INV_MIN = 20_000  # min block size to list an internal inversion
CLEAN_FRAC = 0.90  # contig is "clean" if >= this frac of its aln bp is one chrom + one strand


def chrom_label(name):
    """Short chromosome tag from a messy FASTA header."""
    m = re.search(r"(chr[0-9IVXab]+)", name, re.IGNORECASE)
    if m:
        return m.group(1)
    return name.split("_")[-1] if "_" in name else name


def second_bp_g(d):
    mt = max(d, key=d.get)
    return sum(v for k, v in d.items() if k != mt)


def analyse(paf, min_aln, ref_mb):
    rows = read_paf(paf, min_aln)
    strain = os.path.basename(paf).split("_vs_")[0]
    if not rows:
        return {"strain": strain, "empty": True}

    ref_name = os.path.basename(paf).split("_vs_")[1].replace(".paf", "")
    tlen = {r["t"]: r["tlen"] for r in rows}
    ref_size = ref_mb * 1e6 if ref_mb else sum(tlen.values())

    total_aln = sum(r["te"] - r["ts"] for r in rows)
    ident = sum(r["nmatch"] for r in rows) / sum(r["alen"] for r in rows) * 100

    # reference coverage (merge per chrom)
    byt = collections.defaultdict(list)
    for r in rows:
        byt[r["t"]].append((r["ts"], r["te"]))
    covered = sum(merged_len(v) for v in byt.values())

    # per query contig
    q_rows = collections.defaultdict(list)
    for r in rows:
        q_rows[r["q"]].append(r)

    raw_rev = sum(r["te"] - r["ts"] for r in rows if r["strand"] == "-")

    internal_inv = []  # (chrom, ts, te, span, q)
    transloc = []  # (q, total_bp, {chrom: bp})
    clean_contigs = 0
    for q, rr in q_rows.items():
        bp_by_t = collections.defaultdict(int)
        for r in rr:
            bp_by_t[r["t"]] += r["te"] - r["ts"]
        main_t = max(bp_by_t, key=bp_by_t.get)
        q_tot = sum(bp_by_t.values())
        # dominant strand on main chrom
        strand_bp = collections.defaultdict(int)
        for r in rr:
            if r["t"] == main_t:
                strand_bp[r["strand"]] += r["te"] - r["ts"]
        dom_strand = max(strand_bp, key=strand_bp.get)
        # clean?
        if bp_by_t[main_t] / q_tot >= CLEAN_FRAC and strand_bp[dom_strand] / bp_by_t[main_t] >= CLEAN_FRAC:
            clean_contigs += 1
        # internal inversions: blocks on main chrom against contig's own orientation
        for r in rr:
            if r["t"] == main_t and r["strand"] != dom_strand:
                span = r["te"] - r["ts"]
                if span >= INV_MIN:
                    internal_inv.append((main_t, r["ts"], r["te"], span, q))
        # translocation / chimera: substantial bp to a non-main chrom
        others = {t: bp for t, bp in bp_by_t.items() if t != main_t and bp >= TRANSLOC_MIN}
        if others:
            # query-coord midpoint per ref chrom -> tandem (structural join) vs overlapping (segdup)
            qspan = {}
            for t in list(others) + [main_t]:
                cs = [(r["qs"], r["qe"]) for r in rr if r["t"] == t]
                qspan[t] = (min(s for s, _ in cs), max(e for _, e in cs))
            m0, m1 = qspan[main_t]
            kinds = {}
            for t in others:
                s, e = qspan[t]
                overlap = max(0, min(e, m1) - max(s, m0))
                kinds[t] = "segdup" if overlap > 0.5 * (e - s) else "join"
            transloc.append((q, q_tot, dict([(main_t, bp_by_t[main_t])] + list(others.items())), kinds, qspan))

    internal_inv.sort(key=lambda x: -x[3])
    inv_bp = sum(x[3] for x in internal_inv)

    transloc_bp = sum(second_bp_g(d) for _, _, d, _, _ in transloc)
    # a "join" is a candidate real structural difference; "segdup" is a paralogous/repeat hit
    joins = [
        t
        for t in transloc
        if any(k == "join" for k in t[3].values())
        and sum(v for kk, v in t[2].items() if t[3].get(kk) == "join") >= 100_000
    ]

    # per ref chrom
    chrom_cov = {}
    for t in sorted(tlen, key=lambda t: -tlen[t]):
        cov = merged_len(byt.get(t, []))
        ncontig = len({r["q"] for r in rows if r["t"] == t and r["te"] - r["ts"] >= TRANSLOC_MIN})
        chrom_cov[t] = (cov, tlen[t], ncontig)

    # verdict
    inv_pct = inv_bp / total_aln * 100
    issues = []
    if inv_bp >= 0.005 * total_aln and internal_inv:
        issues.append(f"{len(internal_inv)} inversion block(s) ≥{INV_MIN // 1000}kb ({inv_bp / 1e6:.2f} Mb)")
    if joins:
        issues.append(f"{len(joins)} scaffold(s) joining ≥2 ref chroms end-to-end")
    verdict = "Collinear — 1:1, no rearrangements" if not issues else "; ".join(issues)

    return {
        "strain": strain,
        "ref": ref_name,
        "empty": False,
        "n_aln": len(rows),
        "n_qcontig": len(q_rows),
        "clean": clean_contigs,
        "covered": covered,
        "ref_size": ref_size,
        "cov_pct": covered / ref_size * 100,
        "total_aln": total_aln,
        "ident": ident,
        "raw_rev_pct": raw_rev / total_aln * 100,
        "inv_bp": inv_bp,
        "inv_pct": inv_pct,
        "inv_blocks": internal_inv,
        "transloc": transloc,
        "transloc_bp": transloc_bp,
        "chrom_cov": chrom_cov,
        "verdict": verdict,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("synteny_dir")
    ap.add_argument("--min-aln", type=int, default=2000)
    ap.add_argument(
        "--ref-mb",
        type=float,
        default=None,
        help="reference genome size in Mb for %% coverage (default: sum of aligned chrom lengths)",
    )
    ap.add_argument("--md", default=None)
    ap.add_argument("--tsv", default=None)
    a = ap.parse_args()

    pafs = sorted(glob.glob(os.path.join(a.synteny_dir, "*_vs_*.paf")))
    if not pafs:
        sys.exit(f"no *_vs_*.paf in {a.synteny_dir}")

    res = [analyse(p, a.min_aln, a.ref_mb) for p in pafs]
    res = [r for r in res if not r.get("empty")]
    ref = res[0]["ref"]

    md = []
    md.append(f"### Synteny vs {ref}  (minimap2 -x asm10, alignments ≥{a.min_aln} bp)\n")
    md.append(
        f"| Strain | Query pieces (clean) | {ref} covered | Aln identity¹ | Raw % rev² | Internal inv.³ | Verdict |"
    )
    md.append("| :--- | ---: | ---: | ---: | ---: | ---: | :--- |")
    for r in res:
        md.append(
            f"| **{r['strain']}** | {r['n_qcontig']} ({r['clean']}) "
            f"| {r['covered'] / 1e6:.1f} Mb ({r['cov_pct']:.1f}%) "
            f"| {r['ident']:.1f}% "
            f"| {r['raw_rev_pct']:.0f}% "
            f"| {r['inv_bp'] / 1e6:.2f} Mb ({r['inv_pct']:.1f}%) "
            f"| {r['verdict']} |"
        )
    md.append("")
    md.append(
        "¹ minimap2 gap-compressed identity over aligned regions only; depressed by divergent repeats — whole-genome ANI is the better divergence metric.  "
    )
    md.append(
        "² fraction of aligned bp on the reverse strand *before* per-contig orientation correction; for scaffold assemblies this is ~coin-flip scaffold orientation, not biology.  "
    )
    md.append(
        f"³ alignment blocks ≥{INV_MIN // 1000} kb that run against their own query contig's dominant orientation on the same reference chromosome — the orientation-corrected inversion signal.\n"
    )

    # detail
    md.append("#### Per-strain detail\n")
    for r in res:
        md.append(
            f"**{r['strain']}** — {r['n_aln']} alignments, {r['n_qcontig']} query pieces "
            f"({r['clean']} clean 1:1), {r['covered'] / 1e6:.2f} Mb of {r['ref']} covered "
            f"({r['cov_pct']:.1f}%), identity {r['ident']:.1f}%.  "
        )
        if r["inv_blocks"]:
            top = r["inv_blocks"][:6]
            md.append(
                "  Inversions: "
                + "; ".join(f"{chrom_label(c)}:{s:,}-{e:,} ({sp / 1000:.0f} kb)" for c, s, e, sp, q in top)
                + ".  "
            )
        else:
            md.append("  Inversions: none ≥" + f"{INV_MIN // 1000} kb.  ")
        big_tl = [t for t in r["transloc"] if second_bp_g(t[2]) >= 50_000]
        if big_tl:
            for q, tot, d, kinds, qspan in sorted(big_tl, key=lambda x: -x[1])[:6]:
                mt = max(d, key=d.get)
                parts = ", ".join(
                    f"{chrom_label(k)} {v / 1000:.0f} kb" + ("" if k == mt else f" [{kinds[k]}]")
                    for k, v in sorted(d.items(), key=lambda x: -x[1])
                )
                md.append(f"  Multi-chrom scaffold {q}: {parts}.  ")
        else:
            md.append("  Multi-chromosome scaffolds: none ≥50 kb.  ")
        md.append("")

    out = "\n".join(md)
    print(out)
    if a.md:
        with open(a.md, "w") as fh:
            fh.write(out)
        print(f"\n[wrote {a.md}]", file=sys.stderr)
    if a.tsv:
        with open(a.tsv, "w") as fh:
            fh.write(
                "strain\tref\tn_aln\tn_query_pieces\tclean_pieces\tref_covered_Mb\tref_covered_pct\t"
                "aln_identity_pct\traw_reverse_pct\tinternal_inv_Mb\tinternal_inv_pct\ttransloc_bp\tverdict\n"
            )
            for r in res:
                fh.write(
                    f"{r['strain']}\t{r['ref']}\t{r['n_aln']}\t{r['n_qcontig']}\t{r['clean']}\t"
                    f"{r['covered'] / 1e6:.3f}\t{r['cov_pct']:.2f}\t{r['ident']:.2f}\t{r['raw_rev_pct']:.2f}\t"
                    f"{r['inv_bp'] / 1e6:.3f}\t{r['inv_pct']:.3f}\t{r['transloc_bp']}\t{r['verdict']}\n"
                )
        print(f"[wrote {a.tsv}]", file=sys.stderr)


if __name__ == "__main__":
    main()
