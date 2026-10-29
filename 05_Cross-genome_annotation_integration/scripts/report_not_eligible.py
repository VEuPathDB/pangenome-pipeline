#!/usr/bin/env python3
"""
Report of the genes that could not become panel genes in step 05 (not_eligible.tsv), one row per gene, with the
features behind the decision and an evidence score (nothing is changed).

Per gene:
  outcome            kept_other_model: another member of its basket won and this genome's final GFF has that
                     model, so the gene is kept (the score is about its own model, which was dropped)
                     gone: this genome's final GFF has no model of its basket (no winner, or the winner was removed
                     later or not placed in this genome)
  reason             why step 05 did not let it win (the first rule it failed)
  protein_aa, psauron, has_domain, own_orf (valid / the defect of its own model), valid_orf_genomes (genomes where a
  transfer of it has a valid ORF), transfer_onto_reference
  overlap_gene, overlap_strand, overlap_bp
                     the final-GFF gene of the same genome that it overlaps most (members of its own basket's winners
                     are ignored), i.e. what it would clash with if it had been kept

Evidence score (one point per sign that the gene is not real, 0-4):
  short          protein shorter than --min-aa
  low_psauron    PSAURON below step 05's cut-off for it: --psauron-min-in-ref (0.9) when its sequence lies in the
                 reference (transfer_onto_reference=yes), --psauron-min (0.7) otherwise; or no score
  no_valid_orf   no valid ORF in any genome (own model or transfers)
  overlap        overlaps a gene kept in the final GFF (any strand)
has_domain is shown but not scored (a domain is a sign FOR the gene).

Outputs (--out-dir): not_eligible_report.tsv (sorted by score, highest first), not_eligible_report_summary.tsv and,
for the genes that are gone (really eliminated), per genome:
  <label>_eliminated_reasons.png   the genes by the combination of signs they have (e.g. "PSAURON only")
  <label>_eliminated_psauron.png   histogram of their PSAURON, one panel per genome
  <label>_eliminated_psauron_only.png   the same for the genes eliminated by PSAURON alone

Usage:
    python3 report_not_eligible.py --result-dir basket_annotation --native-dir basket_annotation/input/native \\
        --out-dir basket_annotation/not_eligible_report
"""
import argparse
import csv
import glob
import os
import re
from bisect import bisect_right
from collections import Counter, defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
ATTR_RE = re.compile(r"(?:^|;)([^=;]+)=([^;]*)")


def read_genes(path):
    """genes of a GFF: ID -> (contig, start, end, strand, attributes); a gene is a feature with an ID and no Parent"""
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or line.startswith("#") or "Parent=" in f[8]:
            continue
        m = ID_RE.search(f[8])
        if m:
            out[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6], dict(ATTR_RE.findall(f[8])))
    return out


def genome_of(path, suffix=""):
    return re.sub(suffix + r"\.gff3?$", "", os.path.basename(path))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--result-dir", required=True, help="step 05 output (not_eligible.tsv, baskets_v2.tsv, decide/, "
                                                        "gff_prefer_panel_model_final/)")
    ap.add_argument("--native-dir", required=True, help="the genomes' own GFFs, as step 05 read them")
    ap.add_argument("--final-gff-dir", default=None, help="default: <result-dir>/gff_prefer_panel_model_final")
    ap.add_argument("--min-aa", type=int, default=100)
    ap.add_argument("--psauron-min", type=float, default=0.7, help="step 05's PSAURON cut-off")
    ap.add_argument("--psauron-min-in-ref", type=float, default=0.9,
                    help="step 05's PSAURON cut-off for a gene whose sequence lies in the reference")
    ap.add_argument("--label", default="Toxo")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    final_dir = a.final_gff_dir or os.path.join(a.result_dir, "gff_prefer_panel_model_final")

    ne = list(csv.DictReader(open(os.path.join(a.result_dir, "not_eligible.tsv")), delimiter="\t"))
    winners_of = {}
    for r in csv.DictReader(open(os.path.join(a.result_dir, "baskets_v2.tsv")), delimiter="\t"):
        ws = [w for w in r["winners"].split(";") if w]
        for m in r["members"].split(","):
            winners_of[m] = ws
    own_bad = {}
    p = os.path.join(a.result_dir, "decide", "own_invalid_orf.tsv")
    if os.path.exists(p):
        for r in csv.DictReader(open(p), delimiter="\t"):
            d = [k for k, ok in (("start_lost", r["start_ok"]), ("stop_lost", r["stop_ok"]), ("frameshift", r["frame_ok"]))
                 if ok == "False"]
            if r["internal_stops"] not in ("", "0"):
                d.append(f"internal_stops={r['internal_stops']}")
            own_bad[r["gene"]] = ",".join(d) or "invalid"

    native = {}
    for f in glob.glob(os.path.join(a.native_dir, "*.gff3")) + glob.glob(os.path.join(a.native_dir, "*.gff")):
        native[genome_of(f)] = read_genes(f)
    final = {}   # genome -> contig -> sorted [(start, end, strand, id, basket)]
    kept = set()  # (genome, panel gene) with a model in the final GFF
    for f in glob.glob(os.path.join(final_dir, "*_basket_annotation.gff")):
        by = defaultdict(list)
        for gid, (c, s, e, st, at) in read_genes(f).items():
            by[c].append((s, e, st, gid, at.get("basket", "")))
            kept.add((genome_of(f, "_basket_annotation"), at.get("basket", "")))
        final[genome_of(f, "_basket_annotation")] = {c: sorted(v) for c, v in by.items()}

    rows = []
    for r in ne:
        g, gen = r["gene"], r["genome"]
        ws = winners_of.get(g, [])
        loc = native.get(gen, {}).get(g)
        ov = ("", "", 0)
        if loc and gen in final:
            c, s, e, st, _ = loc
            lst = final[gen].get(c, [])
            i = bisect_right(lst, (e + 1,))
            best = None
            for fs, fe, fst, fid, fb in lst[max(0, i - 200):i]:
                if fe < s or fs > e or fb in ws:
                    continue
                bp = min(e, fe) - max(s, fs) + 1
                if not best or bp > best[2]:
                    best = (fid + (f" ({fb})" if fb else ""), "same" if fst == st else "opposite", bp)
            ov = best or ov
        try:
            aa = int(r["protein_aa"])
        except ValueError:
            aa = None
        ps = float(r["psauron"]) if r["psauron"] not in ("", "-", "NA") else None
        vg = int(r["valid_orf_genomes"]) if r["valid_orf_genomes"].isdigit() else 0
        flags = [k for k, t in (("short", aa is not None and aa < a.min_aa),
                                ("low_psauron", ps is None or ps < (a.psauron_min_in_ref
                                                                     if r["transfer_onto_reference"] == "yes"
                                                                     else a.psauron_min)),
                                ("no_valid_orf", vg == 0),
                                ("overlap", bool(ov[0]))) if t]
        rows.append([len(flags), ",".join(flags), g, gen, "kept_other_model" if any((gen, w) in kept for w in ws) else "gone", ";".join(ws), r["reason"],
                     "" if aa is None else aa, "" if ps is None else f"{ps:.3f}",
                     "yes" if r["has_domain_not_fusion"] == "1" else "no", own_bad.get(g, "valid"), vg,
                     r["transfer_onto_reference"], ov[0], ov[1], ov[2] or "",
                     f"{loc[0]}:{loc[1]}-{loc[2]}({loc[3]})" if loc else ""])

    rows.sort(key=lambda x: (-x[0], x[4] != "gone", x[3], x[2]))
    os.makedirs(a.out_dir, exist_ok=True)
    head = ["score", "flags", "gene", "genome", "outcome", "basket_winners", "reason", "protein_aa", "psauron",
            "has_domain", "own_orf", "valid_orf_genomes", "transfer_onto_reference", "overlap_gene", "overlap_strand",
            "overlap_bp", "location"]
    with open(os.path.join(a.out_dir, "not_eligible_report.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for x in rows:
            fh.write("\t".join(map(str, x)) + "\n")

    c = Counter((x[0], x[4]) for x in rows)
    fl = Counter(f for x in rows for f in x[1].split(",") if f)
    with open(os.path.join(a.out_dir, "not_eligible_report_summary.tsv"), "w") as fh:
        fh.write("score\tgone\tkept_other_model\ttotal\n")
        for s in range(4, -1, -1):
            fh.write(f"{s}\t{c[(s, 'gone')]}\t{c[(s, 'kept_other_model')]}\t"
                     f"{c[(s, 'gone')] + c[(s, 'kept_other_model')]}\n")
        fh.write("\nflag\tgenes\n")
        for k in ("short", "low_psauron", "no_valid_orf", "overlap"):
            fh.write(f"{k}\t{fl[k]}\n")
        fh.write(f"\nwithout_location\t{sum(1 for x in rows if not x[16])}\n")
        fh.write(f"with_domain\t{sum(1 for x in rows if x[9] == 'yes')}\n")
    print(open(os.path.join(a.out_dir, "not_eligible_report_summary.tsv")).read())
    plots([x for x in rows if x[4] == "gone"], a)


FLAG_LABEL = {"low_psauron": "PSAURON", "short": "< 100 aa", "no_valid_orf": "no valid ORF", "overlap": "overlap"}


def plots(gone, a):
    """per genome: the eliminated genes by combination of signs, and the histogram of their PSAURON"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    # every gene with all its signs as recorded: one sign only (PSAURON, < 100 aa or overlap: what alone eliminated
    # the gene), each pair of signs, and "other combinations" = 3 or more signs, no valid ORF as the only sign, and
    # the pairs "< 100 aa + overlap" and "< 100 aa + no valid ORF"
    SIGNS = ("low_psauron", "short", "no_valid_orf", "overlap")
    sig = lambda x: [k for k in SIGNS if k in x[1].split(",")]

    def combo(x):
        fl = [FLAG_LABEL[k] for k in sig(x)]
        return " + ".join(fl) + (" only" if len(fl) == 1 else "")

    MIX = "other combinations"
    to_mix = (["no_valid_orf"], ["short", "overlap"], ["short", "no_valid_orf"])
    cat = lambda x: MIX if len(sig(x)) >= 3 or sig(x) in to_mix else combo(x)
    shown = gone
    genomes = sorted({x[3] for x in gone})
    tot = Counter(cat(x) for x in shown)
    singles = [f"{FLAG_LABEL[k]} only" for k in ("low_psauron", "short", "overlap")]
    pairs = sorted((k for k in tot if " + " in k), key=lambda k: -tot[k])
    cats = [k for k in singles if tot[k]] + pairs + ([MIX] if tot[MIX] else [])
    cnt = Counter((x[3], cat(x)) for x in shown)
    single_col = {"PSAURON only": "#1f5f8b", "< 100 aa only": "#d35400", "overlap only": "#7b4f9d"}
    pair_cols = ["#8fb8d0", "#f0b27a", "#c3a6d6", "#e6a39b", "#a8d5ba", "#d9c27e"]
    colors = [single_col.get(k) or (pair_cols[pairs.index(k) % len(pair_cols)] if k in pairs else "#9a9a9a")
              for k in cats]

    fig, ax = plt.subplots(figsize=(12, 0.42 * len(genomes) + 2.4))
    left = np.zeros(len(genomes))
    for c, col in zip(cats, colors):
        v = np.array([cnt[(g, c)] for g in genomes])
        ax.barh(genomes, v, left=left, color=col, edgecolor="white", lw=0.6, label=f"{c} ({tot[c]:,})")
        for i, (l, n) in enumerate(zip(left, v)):
            if n >= 15:
                ax.text(l + n / 2, i, str(n), ha="center", va="center", fontsize=7,
                        color="white" if c in single_col else "#222")
        left += v
    for i, n in enumerate(left):
        ax.text(n, i, f" {int(n)}", va="center", fontsize=8)
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("genes eliminated (no model of them in the genome's final GFF)")
    leg = ax.legend(title="signs against the gene (total)", frameon=False, fontsize=8, title_fontsize=8.5,
                    loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title(f"{a.label}: genes eliminated by the step 05 quality rules ({len(gone):,})", loc="left", fontsize=12)
    fig.savefig(os.path.join(a.out_dir, f"{a.label}_eliminated_reasons.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    # PSAURON histograms: all eliminated genes, and the ones eliminated by PSAURON alone (same sign set as the
    # "PSAURON only" bar of the reasons figure)
    hist_psauron(gone, genomes, a, "the genes eliminated by the step 05 quality rules",
                 f"{a.label}_eliminated_psauron.png")
    hist_psauron([x for x in gone if cat(x) == "PSAURON only"], genomes, a,
                 "the genes eliminated by PSAURON only", f"{a.label}_eliminated_psauron_only.png")


def hist_psauron(rows, genomes, a, what, name):
    """one PSAURON histogram per genome (10 bins of 0.1), with the counts of the bins from 0.7 up"""
    import matplotlib.pyplot as plt
    import numpy as np
    nc = 4
    nr = -(-len(genomes) // nc)
    fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 2.3 * nr), sharex=True, sharey=True)
    bins = np.linspace(0, 1, 11)   # 10 bins of 0.1
    for axx, g in zip(axs.flat, genomes):
        ps = [float(x[8]) for x in rows if x[3] == g and x[8] != ""]
        n, _, _ = axx.hist(ps, bins=bins, color="#2a6f97", edgecolor="white", lw=0.5)
        for lo, k in zip(bins[:-1], n):   # counts on the bins from 0.7 up (the borderline ones)
            if lo >= 0.7 - 1e-9:
                axx.text(lo + 0.05, k, str(int(k)), ha="center", va="bottom", fontsize=7.5)
        for t in (a.psauron_min, a.psauron_min_in_ref):
            axx.axvline(t, color="#c0392b", ls="--", lw=0.9)
        axx.set_title(f"{g} ({len(ps)})", fontsize=9, loc="left")
        axx.tick_params(labelbottom=True)   # the PSAURON scale on every panel, not only the bottom row
        axx.spines[["top", "right"]].set_visible(False)
    for axx in list(axs.flat)[len(genomes):]:
        axx.set_visible(False)
    fig.supxlabel(f"PSAURON (dashed: step 05 cut-offs {a.psauron_min} and {a.psauron_min_in_ref})", fontsize=10)
    fig.supylabel("genes eliminated", fontsize=10)
    fig.suptitle(f"{a.label}: PSAURON of {what}, by genome ({len(rows):,} genes)", x=0.01, ha="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out_dir, name), dpi=150, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()
