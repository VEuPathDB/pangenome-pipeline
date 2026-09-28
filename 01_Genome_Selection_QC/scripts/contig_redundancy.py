#!/usr/bin/env python3
"""
Report-only contig-redundancy screen for step 01.

Goal: find small contigs that are really duplicated pieces of the chromosomes
(uncollapsed repeats, haplotig-style overlaps, rDNA copies ...). Left in the
assembly they inflate the pangenome with redundant genes.

For every genome, one row per contig, three questions:
  1. Placement  - how much of the contig aligns to the reference (the first genome in
                  genome_order.txt), and to which chromosome  (from  synteny/<label>_vs_<ref>.paf).
  2. Redundancy - how much of the contig is a near-duplicate of a LONGER contig
                  in the same assembly    (from  contig_redundancy/<label>.self.paf,
                  a  minimap2 -xasm5 -DP  self-alignment, the same call purge_dups uses).
  3. Gene load  - how many genes are annotated on it   (from  gff_originals/).

It writes tables + a summary and DROPS NOTHING. Review
contig_redundancy/REPORT.md, then add the contigs you decide to remove to
01_Genome_Selection_QC/drop_contigs.txt by hand and run drop_contigs_slurm.sh
(that tool strips them from genomes/ and gff_originals/ together).

    contig_redundancy.py <project_dir> [options]
"""
import argparse
import collections
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.intervals import merged_len


# ---------------------------------------------------------------- helpers -----
def chrom_label(name):
    """Short chromosome tag from a messy FASTA header (Chr1_A_fumigatus_Af293 -> Chr1)."""
    m = re.search(r"(chr[0-9IVXab]+)", name, re.I)
    if m:
        return m.group(1)
    return name.split("_")[0] if "_" in name else name


def read_fai(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 2:
                d[f[0]] = int(f[1])
    return d


def paf_rows(path, min_aln, min_id):
    """Yield dict rows from a PAF, filtered by block length and gap-compressed identity."""
    out = []
    if not path or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            alen = int(f[10])
            if alen < min_aln:
                continue
            nm = int(f[9])
            if alen and nm / alen < min_id:
                continue
            out.append(
                dict(q=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]),
                     strand=f[4], t=f[5], tlen=int(f[6]),
                     ts=int(f[7]), te=int(f[8]), nm=nm, alen=alen)
            )
    return out


# Gene-level feature types (GFF3 column 3) counted as "a gene" here. NCBI
# GFFs use plain "gene" for everything (biotype is a separate attribute);
# VEuPathDB/FungiDB GFFs (e.g. AF293, A1160, A1163 in this project) commonly
# use "protein_coding_gene" instead, with "pseudogene"/"tRNA_gene"/etc for
# the rest - a bare f[2] == "gene" check recognises none of those and
# silently returns 0 genes for every contig in such a genome. That 0 would
# then read as "safe to drop, gene-free" in the redundancy report, so this
# undercount is not cosmetic - it can make a gene-carrying contig look empty.
GENE_FEATURE_TYPES = {
    "gene", "protein_coding_gene", "pseudogene", "ncRNA_gene",
    "tRNA_gene", "rRNA_gene", "snoRNA_gene", "snRNA_gene", "miRNA_gene",
}


def gene_counts(gff):
    """{seqid: n_gene_features}  or  None if no GFF."""
    if not gff or not os.path.exists(gff):
        return None
    c = collections.Counter()
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 8:
                continue
            if f[2] in GENE_FEATURE_TYPES:
                c[f[0]] += 1
    return c


# ---------------------------------------------------------------- analysis ----
def analyse(label, proj, a):
    fai = os.path.join(proj, a.self_dir, f"{label}.fai")
    if not os.path.exists(fai):
        sys.stderr.write(f"  ! {label}: no {fai} - run the SLURM step first, skipping\n")
        return None
    clen = read_fai(fai)

    # --- 1. placement vs the reference -------------------------------------
    ref_paf = os.path.join(proj, a.synteny_dir, f"{label}_vs_{a.ref}.paf")
    is_ref = not os.path.exists(ref_paf) and label == a.ref
    place = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in paf_rows(ref_paf, a.min_aln, a.ref_min_id):
        place[r["q"]][r["t"]].append((r["qs"], r["qe"]))
    pct_ref, dom_chr = {}, {}
    for c in clen:
        allint = [iv for chrom in place.get(c, {}).values() for iv in chrom]
        pct_ref[c] = 100.0 * merged_len(allint) / clen[c] if clen[c] else 0.0
        if place.get(c):
            bp = {chrom: merged_len(v) for chrom, v in place[c].items()}
            dom_chr[c] = max(bp, key=bp.get)
        else:
            dom_chr[c] = "-"

    # --- 2. self-redundancy ----------------------------------------------
    self_paf = os.path.join(proj, a.self_dir, f"{label}.self.paf")
    dup_int = collections.defaultdict(list)     # contig -> [(s,e)] covered by a longer contig
    covered_by = collections.defaultdict(set)   # contig -> {longer contig names}
    for r in paf_rows(self_paf, a.min_aln, a.dup_min_id):
        q, t = r["q"], r["t"]
        if q == t:
            continue
        lq = clen.get(q, r["qlen"])
        lt = clen.get(t, r["tlen"])
        # attribute the covered stretch to the SHORTER contig (tie -> later name)
        if (lq, q) < (lt, t):
            dup_int[q].append((r["qs"], r["qe"]))
            covered_by[q].add(t)
        else:
            dup_int[t].append((r["ts"], r["te"]))
            covered_by[t].add(q)
    dup_frac = {
        c: (100.0 * merged_len(dup_int[c]) / clen[c] if clen[c] else 0.0)
        for c in clen
    }

    # --- 3. genes -------------------------------------------------------
    gff = os.path.join(proj, a.gff_dir, f"{label}.gff")
    if not os.path.exists(gff):
        gff = os.path.join(proj, a.gff_dir, f"{label}.gff3")
    gc = gene_counts(gff)

    # --- classify -----------------------------------------------------
    big = a.chrom_min_mb * 1e6
    rows = []
    for c in sorted(clen, key=lambda x: -clen[x]):
        L = clen[c]
        ng = "NA" if gc is None else gc.get(c, 0)
        if L >= big:
            cls = "large"
        elif dup_frac[c] >= a.dup_frac:
            cls = "small-redundant"
        elif pct_ref[c] >= a.placed_pct or is_ref:
            cls = "small-placed"
        else:
            cls = "small-unplaced"
        rows.append(dict(
            contig=c, length=L, cls=cls,
            ref_chr=chrom_label(dom_chr[c]) if dom_chr[c] != "-" else "-",
            pct_ref=pct_ref[c], dup_frac=dup_frac[c],
            covered_by=",".join(sorted(covered_by[c])[:3]),
            n_genes=ng,
        ))
    return dict(label=label, is_ref=is_ref, rows=rows, has_genes=gc is not None)


# ---------------------------------------------------------------- output ------
def summarise(res):
    r = res["rows"]
    large = [x for x in r if x["cls"] == "large"]
    small = [x for x in r if x["cls"] != "large"]
    red = [x for x in r if x["cls"] == "small-redundant"]
    unp = [x for x in r if x["cls"] == "small-unplaced"]

    def gsum(rows):
        if not res["has_genes"]:
            return "NA"
        return sum(x["n_genes"] for x in rows)

    return dict(
        label=res["label"],
        n_contig=len(r), n_large=len(large), n_small=len(small),
        small_bp=sum(x["length"] for x in small),
        n_red=len(red), red_bp=sum(x["length"] for x in red), red_genes=gsum(red),
        n_unp=len(unp), unp_bp=sum(x["length"] for x in unp), unp_genes=gsum(unp),
    )


def write_contig_tsv(res, path):
    with open(path, "w") as fh:
        fh.write("contig\tlength\tclass\tref_chr\tpct_aligned_to_ref\t"
                 "pct_duplicated_in_longer_contig\tcovered_by\tn_genes\n")
        for x in res["rows"]:
            fh.write(f"{x['contig']}\t{x['length']}\t{x['cls']}\t{x['ref_chr']}\t"
                     f"{x['pct_ref']:.1f}\t{x['dup_frac']:.1f}\t{x['covered_by']}\t{x['n_genes']}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("project_dir")
    ap.add_argument("--ref", default=None,
                    help="reference label (default: the first genome in <project>/genome_order.txt)")
    ap.add_argument("--order", default="genome_order.txt")
    ap.add_argument("--self-dir", default="contig_redundancy",
                    help="dir with <label>.fai and <label>.self.paf (rel. to project)")
    ap.add_argument("--synteny-dir", default="synteny")
    ap.add_argument("--gff-dir", default="gff_originals")
    ap.add_argument("--min-aln", type=int, default=1000, help="min PAF block length (bp)")
    ap.add_argument("--dup-min-id", type=float, default=0.95,
                    help="min identity for a self-alignment block to count as duplication")
    ap.add_argument("--ref-min-id", type=float, default=0.80,
                    help="min identity for a vs-reference block to count as placed")
    ap.add_argument("--dup-frac", type=float, default=90.0,
                    help="flag a contig as redundant at >= this %% duplicated in a longer contig")
    ap.add_argument("--placed-pct", type=float, default=50.0,
                    help="small contig counts as 'placed' at >= this %% aligned to the reference")
    ap.add_argument("--chrom-min-mb", type=float, default=1.0,
                    help="contigs >= this size are treated as chromosome-scale ('large')")
    ap.add_argument("--md", default=None, help="write the summary report here")
    a = ap.parse_args()

    proj = os.path.abspath(a.project_dir)
    order = os.path.join(proj, a.order)
    if not os.path.exists(order):
        sys.exit(f"no {order}")
    labels = [ln.split()[0] for ln in open(order)
              if ln.strip() and not ln.lstrip().startswith("#")]
    if a.ref is None:
        if not labels:
            sys.exit(f"no genomes in {order}")
        a.ref = labels[0]  # the reference is the first genome, as everywhere else in the pipeline

    out_dir = os.path.join(proj, a.self_dir)
    os.makedirs(out_dir, exist_ok=True)

    results, summaries = [], []
    for lab in labels:
        res = analyse(lab, proj, a)
        if res is None:
            continue
        write_contig_tsv(res, os.path.join(out_dir, f"{lab}.contigs.tsv"))
        results.append(res)
        summaries.append(summarise(res))

    if not results:
        sys.exit("nothing analysed - run contig_redundancy.slurm.sh first")

    md = []
    md.append("### Contig-redundancy screen (report only)\n")
    md.append("Self-alignment: `minimap2 -xasm5 -DP <asm> <asm>` (same call `purge_dups` uses). "
              f"A contig is flagged **small-redundant** when ≥ {a.dup_frac:.0f}% of its length "
              f"aligns (≥ {a.dup_min_id*100:.0f}% identity, ≥ {a.min_aln} bp) to a *longer* "
              "contig of the same assembly. **small-unplaced** = a small contig with "
              f"< {a.placed_pct:.0f}% of its length aligned to {a.ref} and not flagged redundant. "
              "Nothing is removed here — decide per contig from the tables below.\n")
    md.append(f"| Genome | Contigs | large (≥{a.chrom_min_mb:g} Mb) | small | small bp | "
              "redundant | redundant bp | genes on redundant | unplaced | genes on unplaced |")
    md.append("| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for s in summaries:
        md.append(
            f"| {s['label']} | {s['n_contig']} | {s['n_large']} | {s['n_small']} "
            f"| {s['small_bp']/1e6:.2f} Mb | {s['n_red']} | {s['red_bp']/1e3:.0f} kb "
            f"| {s['red_genes']} | {s['n_unp']} | {s['unp_genes']} |"
        )
    md.append("")

    md.append("#### Flagged contigs per genome\n")
    for res in results:
        flagged = [x for x in res["rows"] if x["cls"] in ("small-redundant", "small-unplaced")]
        if not flagged:
            md.append(f"**{res['label']}** — none flagged.\n")
            continue
        md.append(f"**{res['label']}**\n")
        md.append("| contig | length | class | ref chr | % to ref | % dup in longer | covered by | genes |")
        md.append("| :--- | ---: | :--- | :--- | ---: | ---: | :--- | ---: |")
        for x in sorted(flagged, key=lambda y: (y["cls"], -y["length"])):
            md.append(f"| `{x['contig']}` | {x['length']:,} | {x['cls']} | {x['ref_chr']} "
                      f"| {x['pct_ref']:.0f}% | {x['dup_frac']:.0f}% | {x['covered_by'] or '-'} "
                      f"| {x['n_genes']} |")
        md.append("")

    md.append("#### How to read this\n")
    md.append("- **small-redundant + 0 genes** — safe to drop (pure duplicated sequence).")
    md.append("- **small-redundant + genes > 0** — check the genes: BLASTp the proteins on "
              "that contig against the retained contigs; drop only if every gene has a "
              "≥ 98% identity / ≥ 90% coverage hit elsewhere (otherwise the contig carries "
              "real sequence and should stay).")
    md.append("- **small-unplaced** — not redundant and doesn't map to the reference: could be "
              "genuine accessory content, a divergent subtelomere, or contamination. Inspect "
              "(BLAST vs nt, check gene annotations) before deciding.")
    md.append("- **large** contigs are never auto-flagged, no matter how many there are; a "
              "genome with noticeably more large contigs than the reference's chromosome count "
              "means a chromosome got split, not duplicated — that's fragmentation, handled "
              "separately from redundancy.\n")
    md.append("Per-contig tables: `contig_redundancy/<label>.contigs.tsv`.")

    text = "\n".join(md)
    print(text)
    if a.md:
        open(os.path.join(proj, a.md) if not os.path.isabs(a.md) else a.md, "w").write(text)
        sys.stderr.write(f"\n[wrote {a.md}]\n")


if __name__ == "__main__":
    main()
