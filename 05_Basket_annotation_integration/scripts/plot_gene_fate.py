#!/usr/bin/env python3
"""
Two figures from gene_fate.py's table (the fate of every gene of every genome's own annotation):

1. <label>_gene_balance.png: why each genome has fewer (or more) genes in step 05's final GFF than in its own
   annotation. Every own gene and every final gene is counted once; the balance closes exactly
   (own - losses + gains = final), or the script stops.
   Own annotation = the genes step 05 read + the genes removed before step 05 (--removed-genes, --removed-contigs).
   Losses:
     removed with its contig     before step 05, step 01 (the reasons in brackets)
     TE genes removed            before step 05, step 02 (a reason starting with "TE")
     pieces of one gene          kept genes that share one final model (genes - models), and the pieces left out
                                 next to their gene's model (left_out: "piece of a gene that has its model here")
     eliminated by the quality   not eligible (could not win), and left out because the own model failed the check
     rules                       (PSAURON, < 100 aa, invalid ORF)
     removed at the end of       panel gene removed (overlap rule, curation), one place one gene, and left out
     step 05                     because another model is at the same place
   Gains:
     added where the genome had no gene   final models over no gene of the genome's own annotation
     other gains                 one gene with 2+ final models (models - genes), tandem copies added
2. <label>_eliminated_reasons.png: the genes of the red bar ("eliminated by the quality rules"), by the signs against
   them, each checked on its own (not only the first rule the gene failed):
     PSAURON        score below step 05's cut-off: --psauron-min-in-ref (0.9) when the gene's transfer onto the
                    reference works, --psauron-min (0.7) otherwise (not_eligible.tsv), or no score; for a gene left
                    out by the build, below the cut-off of its check (members.tsv)
     < 100 aa       protein shorter than --min-aa
     no valid ORF   no valid ORF in any genome (not_eligible.tsv), or an invalid own ORF (left out by the build)
     overlap        a final gene of another panel gene lies over it (gene_fate's overlapped_by) One sign only,
   each pair, and "other combinations" (3 or more signs, no valid ORF as the only sign, and the pairs "< 100 aa +
   overlap" and "< 100 aa + no valid ORF").

Removed before step 05 (both optional, repeatable, FILE or FILE:REASON):
  --removed-genes   a table with a header: genome (or genome_label), gene (or gene_id), and a reason column (or
                    the REASON given); a biotype column keeps protein_coding only. Or a GFF per genome, with
                    {genome} in the path (e.g. te/{genome}.eliminated_genes.gff:TE genes (step 02)).
  --removed-contigs a table with a header genome, contig, or lines "genome contig" ('#' = comment), with the REASON
                    given; the genes on these contigs are read from --original-gff-dir (the GFFs before the contigs
                    were removed).
The reference genome is left out of the balance (its own genes are the panel genes).

Outputs (--out-dir): <label>_gene_balance.png/.tsv, <label>_eliminated_reasons.png and <label>_eliminated.tsv (one row
per eliminated gene with its signs).

Usage:
    python3 plot_gene_fate.py --gene-fate basket_annotation/gene_fate_prefer_panel_model.tsv \
        --final-dir basket_annotation/gff_prefer_panel_model_final --result-dir basket_annotation \
        --reference ME49 --label Toxo --out-dir basket_annotation/plots \
        --removed-genes 01_Genome_Selection_QC/contig_size_profile/small_contig_genes.tsv:"contig < 10 kb" \
        --removed-contigs 01_Genome_Selection_QC/redundant_contigs_removed.tsv:"redundant contig" \
        --original-gff-dir gff_originals
"""
import argparse
import csv
import glob
import os
import re
from collections import Counter, defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
BASKET_RE = re.compile(r"(?:^|;)basket=([^;\n]+)")
MODEL_RE = re.compile(r"([^;\s]+) \(")
QUALITY = ("PSAURON", "protein shorter", "invalid ORF", "no valid ORF", "not in members.tsv")


def rows(path):
    return list(csv.DictReader(open(path), delimiter="\t")) if path and os.path.exists(path) else []


def split_arg(x):
    """FILE or FILE:REASON (the reason may hold spaces; a path does not hold ':')"""
    f, _, r = x.partition(":")
    return f, r.strip()


def gff_genes(path):
    """protein-coding genes of a GFF (ID and no Parent, with a CDS): ID -> contig"""
    genes, par, cds = {}, {}, set()
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or line.startswith("#"):
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if i and not p:
            genes[i.group(1)] = f[0]
        elif i and p:
            par[i.group(1)] = p.group(1).split(",")[0]
        if f[2] == "CDS" and p:
            t = p.group(1).split(",")[0]
            cds.add(par.get(t, t))
    return {g: c for g, c in genes.items() if g in cds}


def removed_before(a, genomes):
    """genome -> gene -> reason, for the genes removed before step 05"""
    out = defaultdict(dict)
    for x in a.removed_genes:
        f, reason = split_arg(x)
        if "{genome}" in f:
            for g in genomes:
                p = f.replace("{genome}", g)
                if os.path.exists(p):
                    for gid in gff_genes(p):
                        out[g].setdefault(gid, reason or "removed")
            continue
        for r in rows(f):
            g, gid = r.get("genome") or r.get("genome_label"), r.get("gene") or r.get("gene_id")
            if r.get("biotype", "protein_coding") != "protein_coding" or not g or not gid:
                continue
            out[g].setdefault(gid, reason or r.get("reason") or "removed")
    for x in a.removed_contigs:
        f, reason = split_arg(x)
        contigs = defaultdict(set)
        lines = [l for l in open(f) if l.strip() and not l.startswith("#")]
        head = lines[0].rstrip("\n").split("\t")
        if "genome" in head and "contig" in head:
            for r in csv.DictReader(lines, delimiter="\t"):
                contigs[r["genome"]].add(r["contig"])
        else:
            for l in lines:
                t = l.split()
                if len(t) >= 2:
                    contigs[t[0]].add(t[1])
        for g, cs in contigs.items():
            paths = glob.glob(os.path.join(a.original_gff_dir, g + ".gff*")) if a.original_gff_dir else []
            if not paths:
                print(f"WARNING: no GFF of {g} in --original-gff-dir: the genes on its removed contigs are not counted")
                continue
            for gid, c in gff_genes(paths[0]).items():
                if c in cs:
                    out[g].setdefault(gid, reason or "removed contig")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene-fate", required=True, help="gene_fate.py's table, one mode")
    ap.add_argument("--final-dir", required=True, help="the final GFFs of the same mode")
    ap.add_argument("--result-dir", required=True, help="members.tsv and not_eligible.tsv")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--label", default="panel")
    ap.add_argument("--min-aa", type=int, default=100)
    ap.add_argument("--psauron-min", type=float, default=0.7)
    ap.add_argument("--psauron-min-in-ref", type=float, default=0.9)
    ap.add_argument("--removed-genes", action="append", default=[], help="FILE or FILE:REASON (see above)")
    ap.add_argument("--removed-contigs", action="append", default=[], help="FILE or FILE:REASON (see above)")
    ap.add_argument("--original-gff-dir", default="", help="the GFFs before the contigs were removed")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    fate = rows(a.gene_fate)
    genomes = list(dict.fromkeys(r["genome"] for r in fate if r["genome"] != a.reference))
    removed = removed_before(a, genomes)
    contig_reasons = sorted({x for v in removed.values() for x in v.values() if not x.startswith("TE")})

    # the final models of each genome, and which of them are tandem copies
    models, tandem = defaultdict(set), set()
    for g in genomes:
        for line in open(os.path.join(a.final_dir, f"{g}_basket_annotation.gff")):
            f = line.split("\t")
            if len(f) > 8 and f[2] == "gene":
                i = ID_RE.search(f[8]).group(1)
                models[g].add(i)
                if "tandem_array_copy=True" in f[8]:
                    tandem.add(i)

    # signs of the eliminated genes: members.tsv (own model) and not_eligible.tsv (the winner filter)
    mem = {(r["genome"], r["gene"]): r for r in rows(os.path.join(a.result_dir, "members.tsv"))}
    ne = {(r["genome"], r["gene"]): r for r in rows(os.path.join(a.result_dir, "not_eligible.tsv"))}

    def kind(r):
        """the balance category of a gene that has no model of its own basket"""
        f, d = r["fate"], r["detail"]
        if f == "not eligible":
            return "eliminated"
        if f == "left out":
            if d.startswith("piece"):
                return "piece_left_out"
            if d.startswith("same place"):
                return "removed_end"
            if d.startswith(QUALITY):
                return "eliminated"
            return "other"
        if f in ("panel gene removed", "one place, one gene"):
            return "removed_end"
        return "other"

    out, elim = [], []
    for g in genomes:
        fg = [r for r in fate if r["genome"] == g]
        read = {r["gene"] for r in fg}   # a gene that step 05 read was not removed before it
        rm = {k: v for k, v in removed.get(g, {}).items() if k not in read}
        d = Counter()
        # kept genes and their final models: genes sharing models = pieces joined; a gene with 2+ models = split
        parent = {}

        def find(x):
            while parent.setdefault(x, x) != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        used_models = set()
        for r in fg:
            if r["fate"] != "kept":
                k = kind(r)
                d[k] += 1
                if k == "eliminated":
                    elim.append(r)
                continue
            for m in MODEL_RE.findall(r["detail"]):
                used_models.add(m)
                parent[find("g:" + r["gene"])] = find("m:" + m)
        comp = defaultdict(lambda: [0, 0])
        for x in list(parent):
            comp[find(x)][0 if x.startswith("g:") else 1] += 1
        for n, m in comp.values():
            if n > m:
                d["pieces_joined"] += n - m
            elif m > n:
                d["split"] += m - n
        free = models[g] - used_models
        d["added"] = sum(1 for m in free if m not in tandem)
        d["extra_copies"] = sum(1 for m in free if m in tandem)
        for x in contig_reasons:
            d[f"removed: {x}"] = sum(1 for y in rm.values() if y == x)
        d["removed_te"] = sum(1 for y in rm.values() if y.startswith("TE"))
        d["removed_with_contig"] = len(rm) - d["removed_te"]
        own, final = len(fg) + len(rm), len(models[g])
        losses = (d["removed_with_contig"] + d["removed_te"] + d["pieces_joined"] + d["piece_left_out"]
                  + d["eliminated"] + d["removed_end"] + d["other"])
        gains = d["added"] + d["split"] + d["extra_copies"]
        if own - losses + gains != final:
            raise SystemExit(f"ERROR: {g}: the balance does not close: {own} - {losses} + {gains} != {final}")
        if d["other"]:
            print(f"WARNING: {g}: {d['other']} genes with a fate this script does not know (counted as 'other')")
        out.append(dict(d, genome=g, own=own, final=final))

    keys = (["genome", "own", "final", "removed_with_contig"] + [f"removed: {x}" for x in contig_reasons]
            + ["removed_te", "pieces_joined", "piece_left_out", "eliminated", "removed_end", "other", "added",
               "split", "extra_copies"])
    with open(os.path.join(a.out_dir, f"{a.label}_gene_balance.tsv"), "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for d in out:
            fh.write("\t".join(str(d.get(k, 0)) for k in keys) + "\n")

    # 1. the balance
    contig_label = "removed with its contig before step 05" + (f" ({', '.join(contig_reasons)})" if contig_reasons else "")
    fig_losses = [(contig_label, ["removed_with_contig"], "#5d4037", True),
                  ("TE genes removed (step 02)", ["removed_te"], "#b08968", True),
                  ("pieces of one gene replaced by one model", ["pieces_joined", "piece_left_out"], "#1f5f8b", False),
                  ("eliminated by the step 05 quality rules", ["eliminated"], "#c0392b", False),
                  ("removed at the end of step 05 (overlap rule, curation, same place)", ["removed_end"], "#e6a39b",
                   False),
                  ("other", ["other"], "#bdbdbd", False)]
    fig_gains = [("added where the genome had no gene", ["added"], "#2e8b57"),
                 ("other gains (one gene split in two, tandem copies)", ["split", "extra_copies"], "#a8d5ba")]
    fig, ax = plt.subplots(figsize=(13, 0.5 * len(out) + 2.6))
    y = np.arange(len(out))[::-1]
    left = np.zeros(len(out))
    for lab, ks, col, keep in fig_losses:
        v = np.array([sum(d.get(k, 0) for k in ks) for d in out])
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
        v = np.array([sum(d.get(k, 0) for k in ks) for d in out])
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
    plt.close(fig)

    # 2. the eliminated genes by their signs
    LABEL = {"low_psauron": "PSAURON", "short": f"< {a.min_aa} aa", "no_valid_orf": "no valid ORF",
             "overlap": "overlap"}
    erows = []
    for r in elim:
        k = (r["genome"], r["gene"])
        m, n = mem.get(k, {}), ne.get(k, {})
        why = n.get("reason") or r["detail"]
        signs = []
        if n:   # could not win: the winner filter's values
            aa, ps = int(n.get("protein_aa") or 0), n.get("psauron", "")
            thr = a.psauron_min_in_ref if n.get("transfer_onto_reference") == "yes" else a.psauron_min
            low = ps == "" or float(ps) < thr
            no_orf = n.get("valid_orf_genomes") == "0"
        else:   # left out by the build: the check of its own model (members.tsv)
            aa, ps = int(m.get("protein_aa") or 0), m.get("psauron", "")
            t = re.search(r"PSAURON [\d.]+ < ([\d.]+)", m.get("keep_check", ""))
            thr = float(t.group(1)) if t else a.psauron_min
            low = ps != "" and float(ps) < thr
            no_orf = m.get("own_orf") == "invalid"
        if low:
            signs.append("low_psauron")
        if aa < a.min_aa:
            signs.append("short")
        if no_orf:
            signs.append("no_valid_orf")
        if r.get("overlapped_by", "-") not in ("", "-"):
            signs.append("overlap")
        erows.append((r, signs, aa, ps, why))
    with open(os.path.join(a.out_dir, f"{a.label}_eliminated.tsv"), "w") as fh:
        fh.write("genome\tgene\tlocation\tfate\treason\tsigns\tprotein_aa\tpsauron\toverlapped_by\n")
        for r, signs, aa, ps, why in erows:
            fh.write(f"{r['genome']}\t{r['gene']}\t{r['location']}\t{r['fate']}\t{why}\t{','.join(signs) or '-'}\t"
                     f"{aa}\t{ps}\t{r.get('overlapped_by', '-')}\n")

    SIGNS = ("low_psauron", "short", "no_valid_orf", "overlap")
    MIX = "other combinations"
    to_mix = (["no_valid_orf"], ["short", "overlap"], ["short", "no_valid_orf"])

    def cat(signs):
        s = [k for k in SIGNS if k in signs]
        if len(s) >= 3 or s in to_mix or not s:
            return MIX
        return " + ".join(LABEL[k] for k in s) + (" only" if len(s) == 1 else "")
    tot = Counter(cat(x[1]) for x in erows)
    cnt = Counter((x[0]["genome"], cat(x[1])) for x in erows)
    singles = [f"{LABEL[k]} only" for k in ("low_psauron", "short", "overlap")]
    pairs = sorted((k for k in tot if " + " in k), key=lambda k: -tot[k])
    cats = [k for k in singles if tot[k]] + pairs + ([MIX] if tot[MIX] else [])
    single_col = {"PSAURON only": "#1f5f8b", f"< {a.min_aa} aa only": "#d35400", "overlap only": "#7b4f9d"}
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
    ax.set_xlabel("genes eliminated (no model of their own panel gene in the genome's final GFF)")
    ax.legend(title="signs against the gene (total)", frameon=False, fontsize=8, title_fontsize=8.5,
              loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.set_title(f"{a.label}: genes eliminated by the step 05 quality rules ({len(erows):,})", loc="left", fontsize=12)
    fig.savefig(os.path.join(a.out_dir, f"{a.label}_eliminated_reasons.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    for x in ("gene_balance.png", "gene_balance.tsv", "eliminated_reasons.png", "eliminated.tsv"):
        print(os.path.join(a.out_dir, f"{a.label}_{x}"))


if __name__ == "__main__":
    main()
