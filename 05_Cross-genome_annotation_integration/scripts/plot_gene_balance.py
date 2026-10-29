#!/usr/bin/env python3
"""
Why a genome has fewer (or more) genes after step 05 than in its own annotation: every gene of the own annotation and
every gene of the final GFF (prefer_panel_model) accounted for, per genome. The balance closes exactly
(own - losses + gains = final); the script stops if it does not.

"Own annotation" is the genome's original annotation: the GFF step 05 read plus the genes removed with their contigs
before step 05 (--removed-with-contigs: a table genome, gene, contig, reason, e.g. contigs < 10 kb, apicoplast or
redundant contigs, or TE genes). Two categories, always shown (also when 0): TE genes removed (step 02; a reason
starting with "TE") and genes removed with their contig (every other reason); the table has one column per reason.

Losses (own genes that do not become one final gene each):
  <reason>             removed with its contig before step 05 (one category per reason of --removed-with-contigs)
  eliminated           the own gene failed step 05's rules and has no model in the final GFF: the same genes as the
                       not_eligible report's "gone" (so both figures give the same number)
  pieces joined        2+ own genes (pieces of one gene) became one model of the final GFF: counted as pieces - models
  piece left out       the own gene is in a basket that has a model in this genome, and it overlaps that model or lies
                       next to it (< --near bp): a piece of a gene that is still there, the model does not use it
  copy elsewhere       the same, but farther than --near bp or on another contig: maybe a second copy, no model
  model removed later  its basket's model is not in the final GFF of this genome: step 05 removed it later (panel gene
                       removed from every genome by the overlap rule or curation, or another model took its place)
Only the models in the final GFF count: a model step 05 placed and removed later is neither a gain nor a loss.
  not in step 05 input the own gene is in no basket and not in --removed-with-contigs (should be 0)
Every own gene is put in exactly one category (nothing is computed as a remainder).
Gains (final genes that are not one own gene each):
  added (no own gene)  a model placed where the genome had no gene of its own
  own gene split       one own gene lies under 2+ panel genes: counted as models - own genes
  extra copies         extra tandem array copies

The reference genome is left out (its own genes are the panel genes).

Outputs (--out-dir): <label>_gene_balance.png, <label>_gene_balance.tsv and <label>_copies_without_model.tsv (the
"copy elsewhere" genes, one row each: the gene, where it is, its panel gene, and the model of that panel gene in the
same genome, i.e. the copy that kept the model)

Usage:
    python3 plot_gene_balance.py --result-dir basket_annotation --report-dir basket_annotation/not_eligible_report \\
        --native-dir basket_annotation/input/native \\
        --removed-with-contigs 01_Genome_Selection_QC/genes_removed_with_contigs.tsv \\
        --reference ME49 --label Toxo --out-dir basket_annotation/plots
"""
import argparse
import csv
import glob
import os
import re
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys_path = os.path.dirname(os.path.abspath(__file__))
SPLIT = re.compile(r"[,;]")
LOSSES = [("pieces_joined", "pieces of one gene joined into one model", "#1f5f8b"),
          ("eliminated", "eliminated (step 05 rules)", "#c0392b"),
          ("model_removed_later", "model removed later in step 05 (overlap rule, curation, same place)", "#e6a39b"),
          ("piece_left_out", "piece of a gene left out of its model (gene still there)", "#8fb8d0"),
          ("copy_elsewhere", "copy elsewhere, no model", "#9a9a9a"),
          ("not_in_input", "not in step 05's input", "#4d4d4d")]
GAINS = [("added", "added where the genome had no gene", "#2e8b57"),
         ("split", "one own gene split into 2+ genes", "#a8d5ba"),
         ("extra_copies", "extra tandem copies", "#d9c27e")]


def own_genes(native_dir, genome):
    """protein-coding genes of the own annotation: ID -> (contig, start, end) of genes (ID, no Parent) with a CDS"""
    paths = [p for p in glob.glob(os.path.join(native_dir, genome + ".gff*"))]
    genes, par, cds = {}, {}, set()
    for line in open(paths[0]):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or line.startswith("#"):
            continue
        i, p = re.search(r"(?:^|;)ID=([^;]+)", f[8]), re.search(r"(?:^|;)Parent=([^;]+)", f[8])
        if i and not p:
            dm = re.search(r"(?:^|;)description=([^;]+)", f[8])
            genes[i.group(1)] = (f[0], int(f[3]), int(f[4]), dm.group(1) if dm else "")
        elif i and p:
            par[i.group(1)] = p.group(1).split(",")[0]
        if f[2] == "CDS" and p:
            t = p.group(1).split(",")[0]
            cds.add(par.get(t, t))
    return {g: v for g, v in genes.items() if g in cds}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--result-dir", required=True, help="step 05 output (decide/decisions.tsv, final GFFs)")
    ap.add_argument("--report-dir", required=True, help="report_not_eligible.py output (not_eligible_report.tsv)")
    ap.add_argument("--native-dir", required=True, help="the genomes' own GFFs, as step 05 read them")
    ap.add_argument("--removed-with-contigs", default=None,
                    help="genes removed with their contigs before step 05 (genome, gene, contig, reason)")
    ap.add_argument("--contig-reasons", default=None,
                    help="the reasons of the contig removal, for the legend (default: those of --removed-with-contigs)")
    ap.add_argument("--near", type=int, default=2000, help="piece left out: within this distance of its model")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--label", default="Toxo")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()

    dec = list(csv.DictReader(open(os.path.join(a.result_dir, "decide", "decisions.tsv")), delimiter="\t"))
    rep = list(csv.DictReader(open(os.path.join(a.report_dir, "not_eligible_report.tsv")), delimiter="\t"))
    final_dir = os.path.join(a.result_dir, "gff_prefer_panel_model_final")
    genomes = sorted({r["genome"] for r in dec} - {a.reference})
    winners_of = {}
    for r in csv.DictReader(open(os.path.join(a.result_dir, "baskets_v2.tsv")), delimiter="\t"):
        for m in r["members"].split(","):
            winners_of[m] = [w for w in r["winners"].split(";") if w]
    in_basket = set(winners_of)

    removed = {}   # genome -> gene -> reason
    if a.removed_with_contigs:
        for r in csv.DictReader(open(a.removed_with_contigs), delimiter="\t"):
            removed.setdefault(r["genome"], {})[r["gene"]] = r["reason"]
    reasons = sorted({x for v in removed.values() for x in v.values()},
                     key=lambda x: -sum(1 for v in removed.values() for y in v.values() if y == x))
    contig_reasons = ([x.strip() for x in a.contig_reasons.split(",")] if a.contig_reasons
                      else [x for x in reasons if not x.startswith("TE")])
    losses = [("removed_with_contig", "removed with its contig before step 05"
               + (f" ({', '.join(contig_reasons)})" if contig_reasons else ""), "#5d4037"),
              ("removed_te", "TE genes removed (step 02)", "#b08968")] + LOSSES
    always = {"removed_with_contig", "removed_te"}

    out, copies = [], []
    for g in genomes:
        fb, model_at = Counter(), {}
        for line in open(os.path.join(final_dir, f"{g}_basket_annotation.gff")):
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                b = re.search(r"basket=([^;\n]+)", f[8]).group(1)
                fb[b] += 1
                om = re.search(r"(?:^|;)original_ID=([^;\n]+)", f[8])
                model_at.setdefault(b, []).append((f[0], int(f[3]), int(f[4]), om.group(1) if om else "-"))
        final = sum(fb.values())
        rm = removed.get(g, {})
        loc = own_genes(a.native_dir, g)
        owns = set(loc) - set(rm)
        own = len(owns) + len(rm)
        gone = {r["gene"] for r in rep if r["genome"] == g and r["outcome"] == "gone"} & owns
        # the models of the final GFF and the own genes each one stands for (eliminated genes excluded)
        rows = [r for r in dec if r["genome"] == g and r["decision_prefer_panel_model"] != "absent"
                and fb[r["panel_gene"]] > 0]
        own_of = {r["panel_gene"]: [x for x in SPLIT.split(r["own_genes"]) if x and x not in gone] for r in rows}
        cnt = Counter(x for v in own_of.values() for x in v)
        multi = [v for v in own_of.values() if len(v) > 1]
        in_multi = {x for v in multi for x in v}
        shared = {x for x, n in cnt.items() if n > 1} - in_multi
        one = sum(1 for v in own_of.values() if len(v) == 1 and cnt[v[0]] == 1)
        # models made from one own gene that also makes other models (alone, or as a piece of a joined model)
        split_rows = sum(1 for v in own_of.values() if len(v) == 1 and cnt[v[0]] > 1)
        placed = set(cnt)
        if set(rm) & placed:
            raise SystemExit(f"ERROR: {g}: genes listed as removed before step 05 are used by step 05")
        d = dict(genome=g, own=own, final=final, one_to_one=one, pieces=len(in_multi), models_from_pieces=len(multi),
                 pieces_joined=len(in_multi) - len(multi),
                 eliminated=len(gone),
                 added=sum(1 for v in own_of.values() if not v),
                 split=split_rows - len(shared),
                 extra_copies=sum(n - 1 for p, n in fb.items() if n > 1 and p in own_of)
                 + sum(n for p, n in fb.items() if p not in own_of))
        for x in reasons:
            d[f"removed: {x}"] = sum(1 for y in rm.values() if y == x)
        d["removed_te"] = sum(1 for y in rm.values() if y.startswith("TE"))
        d["removed_with_contig"] = len(rm) - d["removed_te"]
        # own genes in a basket with no own model of theirs: where is their basket's model in this genome?
        where = Counter()
        for x in (owns & in_basket) - placed - gone:
            ms = [m for w in winners_of[x] for m in model_at.get(w, [])]
            c, s0, e0, desc = loc[x]
            near = [m for m in ms if m[0] == c and max(m[1], s0) - min(m[2], e0) < a.near]
            k = "model_removed_later" if not ms else "piece_left_out" if near else "copy_elsewhere"
            where[k] += 1
            if k == "copy_elsewhere":
                same = [m for m in ms if m[0] == c]
                dist = min(max(m[1], s0) - min(m[2], e0) for m in same) if same else ""
                copies.append([g, x, f"{c}:{s0}-{e0}", desc, ";".join(winners_of[x]),
                               ";".join(f"{m[0]}:{m[1]}-{m[2]}" for m in ms), ";".join(m[3] for m in ms),
                               "same contig" if same else "other contig", dist])
        d.update({k: where[k] for k in ("piece_left_out", "copy_elsewhere", "model_removed_later")})
        d["not_in_input"] = len(owns - in_basket - placed - gone)
        extra_own = placed - owns   # own genes named by step 05 but not in this GFF
        if extra_own:
            raise SystemExit(f"ERROR: {g}: {len(extra_own)} own genes of step 05 are not in {a.native_dir}")
        calc = own - sum(d[k] for k, _, _ in losses) + sum(d[k] for k, _, _ in GAINS)
        if calc != final:
            raise SystemExit(f"ERROR: {g}: balance does not close (own {own} -> {calc}, final GFF {final})")
        out.append(d)

    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, f"{a.label}_copies_without_model.tsv"), "w") as fh:
        fh.write("genome\tgene\tlocation\tdescription\tpanel_gene\tmodel_location\tmodel_original_ID\t"
                 "model_contig\tdistance_bp\n")
        for r in copies:
            fh.write("\t".join(map(str, r)) + "\n")
    keys = ["genome", "own", "final", "one_to_one", "pieces", "models_from_pieces"] + \
           [k for k, _, _ in losses] + [f"removed: {x}" for x in reasons] + \
           [k for k, _, _ in GAINS]
    with open(os.path.join(a.out_dir, f"{a.label}_gene_balance.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for d in out:
            fh.write("\t".join(str(d[k]) for k in keys) + "\n")

    # the figure groups the table's columns into few bars (the table keeps every column)
    fig_losses = [(losses[0][1], ["removed_with_contig"], "#5d4037", True),   # with the reasons in brackets
                  ("TE genes removed (step 02)", ["removed_te"], "#b08968", True),
                  ("pieces of one gene replaced by one model", ["pieces_joined", "piece_left_out"], "#1f5f8b", False),
                  ("eliminated by the step 05 quality rules", ["eliminated"], "#c0392b", False),
                  ("removed at the end of step 05 (overlap rule, curation, same place)", ["model_removed_later"],
                   "#e6a39b", False),
                  ("copy elsewhere, no model", ["copy_elsewhere", "not_in_input"], "#bdbdbd", False)]
    fig_gains = [("added where the genome had no gene", ["added"], "#2e8b57"),
                 ("other gains (one gene split in two, tandem copies)", ["split", "extra_copies"], "#a8d5ba")]
    fig, ax = plt.subplots(figsize=(13, 0.5 * len(out) + 2.6))
    y = np.arange(len(out))[::-1]
    left = np.zeros(len(out))
    for lab, ks, col, keep in fig_losses:
        v = np.array([sum(d[k] for k in ks) for d in out])
        if not v.sum() and not keep:
            continue
        ax.barh(y, -v, left=-left, color=col, edgecolor="white", lw=0.6, label=f"- {lab}")
        for yy, l, n in zip(y, left, v):
            if n >= 120:
                ax.text(-(l + n / 2), yy, f"{n:,}", ha="center", va="center", fontsize=7.5,
                        color="white" if col in ("#1f5f8b", "#c0392b", "#5d4037") else "#222")
        left += v
    right = np.zeros(len(out))
    for lab, ks, col in fig_gains:
        v = np.array([sum(d[k] for k in ks) for d in out])
        ax.barh(y, v, left=right, color=col, edgecolor="white", lw=0.6, label=f"+ {lab}")
        for yy, l, n in zip(y, right, v):
            if n >= 120:
                ax.text(l + n / 2, yy, f"{n:,}", ha="center", va="center", fontsize=7.5, color="white")
        right += v
    xmax = right.max() * 1.05
    for yy, d in zip(y, out):
        ax.text(xmax + 40, yy, f"{d['own']:,} -> {d['final']:,}  ({d['final'] - d['own']:+,})", va="center",
                fontsize=8.5)
    ax.axvline(0, color="#222", lw=0.8)
    ax.set_yticks(y, [d["genome"] for d in out])
    ax.set_xlim(-left.max() * 1.05, xmax + 1200)
    ax.set_xlabel("genes lost (left) and gained (right) between the genome's own annotation and the final GFF")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1))
    ax.set_title(f"{a.label}: from each genome's own annotation to step 05's final GFF", loc="left", fontsize=12)
    fig.savefig(os.path.join(a.out_dir, f"{a.label}_gene_balance.png"), dpi=150, bbox_inches="tight")
    print(os.path.join(a.out_dir, f"{a.label}_gene_balance.png"))


if __name__ == "__main__":
    main()
