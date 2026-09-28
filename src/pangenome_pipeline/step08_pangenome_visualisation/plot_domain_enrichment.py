#!/usr/bin/env python3
"""InterPro term enrichment per pangene occupancy category (dot plot).

For each occupancy category (Core / Soft-core / Shell / Cloud, same definition as plot_occupancy.py) and
each InterPro entry, a one-sided Fisher exact test asks whether clusters of that category carry the entry
more often than all the other clusters. p-values are Benjamini-Hochberg adjusted over all tests.
Dot plot: x = category, y = InterPro entry, dot size = gene ratio (clusters in the category with the entry /
clusters in the category with any InterPro entry), dot colour = adjusted p-value.

A cluster carries an InterPro entry if any gene of the cluster (any genome) has it in the Dbxref= of its
CDS lines. Only clusters with at least one InterPro entry are counted (the test universe).

Entries carried by fewer than --min-total clusters of the whole universe (default 5) are not tested, so they do not
count in the Benjamini-Hochberg correction. The cut is on the total size of the entry, which does not depend on how
the entry is split between the categories. --min-count (default 3) is a different filter: it only decides which
significant entries are drawn, and needs at least that many clusters of the category to carry the entry.

Inputs:
  --genes-matrix  pangene_matrix_genes.tr.tab from get_pangenes (rows = clusters, columns = genomes,
                  cells = gene IDs, "-" = absent)
  --gff-dir       GFFs whose gene IDs match the matrix (the ones given to get_pangenes), any *.gff / *.gff3
  --entry-list    InterPro entry.list (IPR id -> name); downloaded if not given
Writes <out-prefix>.png/.svg and <out-prefix>.tsv (all tests).

Needs numpy, scipy, matplotlib.
"""

import argparse
import csv
import glob
import os
import re
import sys
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import false_discovery_control, fisher_exact

from pangenome_pipeline.step08_pangenome_visualisation.plot_occupancy import ORDER, categorize

ENTRY_LIST_URL = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"
ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
DBXREF_RE = re.compile(r"(?:^|;)Dbxref=([^;]+)")


def gene_domains(gff_dir):
    """gene ID -> set of IPR ids, from CDS Dbxref (CDS -> mRNA -> gene)."""
    paths = sorted(glob.glob(os.path.join(gff_dir, "*.gff")) + glob.glob(os.path.join(gff_dir, "*.gff3")))
    if not paths:
        sys.exit(f"ERROR: no .gff/.gff3 files in {gff_dir}")
    out = {}
    for path in paths:
        mrna_gene, mrna_ipr = {}, {}
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                if f[2] == "mRNA" or f[2] == "transcript":
                    i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
                    if i and p:
                        mrna_gene[i.group(1)] = p.group(1)
                elif f[2] == "CDS":
                    p, d = PARENT_RE.search(f[8]), DBXREF_RE.search(f[8])
                    if p and d:
                        ipr = {t.split(":", 1)[1] for t in d.group(1).split(",") if t.startswith("InterPro:")}
                        if ipr:
                            mrna_ipr.setdefault(p.group(1), set()).update(ipr)
        for mrna, ipr in mrna_ipr.items():
            out.setdefault(mrna_gene.get(mrna, mrna), set()).update(ipr)
    return out


def read_clusters(path):
    """cluster -> (occupancy, [gene IDs])."""
    clusters = {}
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    keep = [i for i, g in enumerate(rows[0]) if i > 0 and g.strip()]
    for r in rows[1:]:
        if r[0].startswith("chr"):
            continue
        cells = [r[i] for i in keep if i < len(r) and r[i] not in ("", "-", "NA")]
        if cells:
            clusters[r[0]] = (len(cells), [g for c in cells for g in c.split(",")])
    return clusters, len(keep)


def load_entry_names(path, out_prefix):
    if not path:
        path = out_prefix + "_entry.list"
        if not os.path.exists(path):
            try:
                print(f"Downloading {ENTRY_LIST_URL}", file=sys.stderr)
                urllib.request.urlretrieve(ENTRY_LIST_URL, path)
            except Exception as e:  # noqa: BLE001 - any download failure falls back to IPR ids
                print(f"WARNING: could not download entry.list ({e}); using IPR ids as labels", file=sys.stderr)
                return {}
    names = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0].startswith("IPR"):
                names[f[0]] = f[2]
    return names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--entry-list", default="")
    ap.add_argument("-o", "--out-prefix", required=True)
    ap.add_argument(
        "--core-missing", type=int, default=0, help="core = missing from at most this many genomes (default 0)"
    )
    ap.add_argument(
        "--softcore-missing", type=int, default=1, help="soft-core = missing from at most this many genomes (default 1)"
    )
    ap.add_argument("--cloud-max", type=int, default=2, help="cloud = present in at most this many genomes (default 2)")
    ap.add_argument("--alpha", type=float, default=0.05, help="adjusted p-value cutoff (default 0.05)")
    ap.add_argument(
        "--min-total",
        type=int,
        default=5,
        help="entries carried by fewer clusters than this in the whole universe are not tested (default 5)",
    )
    ap.add_argument(
        "--min-count",
        type=int,
        default=3,
        help="only entries carried by at least this many clusters of the category are drawn (default 3)",
    )
    ap.add_argument("--top", type=int, default=15, help="max entries shown per category (default 15)")
    ap.add_argument(
        "--categories", default="Core,Cloud", help="comma-separated categories to test and plot (default Core,Cloud)"
    )
    ap.add_argument("--title", default="")
    args = ap.parse_args()

    wanted = [c.strip() for c in args.categories.split(",")]
    bad = [c for c in wanted if c not in ORDER]
    if bad:
        sys.exit(f"ERROR: unknown categories {bad}; choose from {ORDER}")

    wanted.sort(key=ORDER.index, reverse=True)  # Cloud first (left), Core last (right), like the occupancy plot
    clusters, n_genomes = read_clusters(args.genes_matrix)
    domains = gene_domains(args.gff_dir)
    print(f"{len(clusters)} clusters, {n_genomes} genomes, {len(domains)} genes with InterPro entries", file=sys.stderr)

    cluster_cat, cluster_ipr = {}, {}
    for cid, (occ, genes) in clusters.items():
        ipr = set().union(*(domains.get(g, set()) for g in genes)) if genes else set()
        if ipr:
            cluster_cat[cid] = categorize(occ, n_genomes, args.core_missing, args.softcore_missing, args.cloud_max)
            cluster_ipr[cid] = ipr
    N = len(cluster_ipr)
    if N == 0:
        sys.exit("ERROR: no cluster has InterPro entries; check that --gff-dir IDs match the matrix gene IDs")
    print(f"{N} clusters with at least one InterPro entry (test universe)", file=sys.stderr)

    total = {}
    per_cat = {c: {} for c in wanted}
    n_cat = {c: 0 for c in wanted}
    for cid, ipr in cluster_ipr.items():
        cat = cluster_cat[cid]
        for t in ipr:
            total[t] = total.get(t, 0) + 1
        if cat in per_cat:
            n_cat[cat] += 1
            for t in ipr:
                per_cat[cat][t] = per_cat[cat].get(t, 0) + 1

    rows = []
    n_tested = sum(1 for K in total.values() if K >= args.min_total)
    print(
        f"{len(total)} InterPro entries in the universe; {n_tested} carried by at least {args.min_total} clusters are tested",
        file=sys.stderr,
    )
    for cat in wanted:
        n = n_cat[cat]
        for t, k in per_cat[cat].items():
            K = total[t]
            if K < args.min_total:
                continue
            _, p = fisher_exact([[k, n - k], [K - k, N - n - (K - k)]], alternative="greater")
            rows.append([cat, t, k, n, k / n, K, p])
    if not rows:
        sys.exit("ERROR: nothing to test (no entry reaches --min-total in the requested categories)")
    padj = false_discovery_control([r[6] for r in rows], method="bh")
    for r, q in zip(rows, padj):
        r.append(q)

    names = load_entry_names(args.entry_list, args.out_prefix)
    with open(args.out_prefix + ".tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            [
                "category",
                "interpro",
                "name",
                "clusters_with_entry_in_category",
                "clusters_in_category",
                "gene_ratio",
                "clusters_with_entry_total",
                "p",
                "p_adjusted",
            ]
        )
        for r in sorted(rows, key=lambda r: (r[0], r[7])):
            w.writerow([r[0], r[1], names.get(r[1], ""), r[2], r[3], f"{r[4]:.5f}", r[5], f"{r[6]:.3g}", f"{r[7]:.3g}"])

    shown = []
    for cat in wanted:
        sig = [r for r in rows if r[0] == cat and r[7] < args.alpha and r[2] >= args.min_count]
        sig.sort(key=lambda r: r[7])
        shown += sig[: args.top]
    if not shown:
        sys.exit(
            f"No entry passes adjusted p < {args.alpha} and >= {args.min_count} clusters; see {args.out_prefix}.tsv"
        )

    labels = {}
    for r in shown:
        labels[r[1]] = names.get(r[1], r[1])
    terms = []
    for cat in wanted:
        terms += [r[1] for r in sorted((r for r in shown if r[0] == cat), key=lambda r: r[4])]
    seen = list(dict.fromkeys(terms))
    ypos = {t: i for i, t in enumerate(reversed(seen))}
    cats_plotted = [c for c in wanted if any(r[0] == c for r in shown)]
    xpos = {c: i for i, c in enumerate(cats_plotted)}

    ratios = np.array([r[4] for r in shown])
    sizes = 30 + 250 * (ratios / ratios.max())
    fig, ax = plt.subplots(figsize=(8.5, 0.28 * len(seen) + 2.5))
    sc = ax.scatter(
        [xpos[r[0]] for r in shown],
        [ypos[r[1]] for r in shown],
        s=sizes,
        c=[r[7] for r in shown],
        cmap="turbo",
        edgecolors="none",
    )
    ax.set_yticks(list(ypos.values()))
    ax.set_yticklabels([labels[t] for t in ypos])
    ax.set_xticks(list(xpos.values()))
    ax.set_xticklabels(cats_plotted)
    ax.set_xlim(-0.6, len(cats_plotted) - 0.4)
    ax.set_ylim(-0.7, len(seen) - 0.3)
    ax.set_xlabel("Pangene occupancy", fontsize=12)
    ax.set_ylabel("InterPro terms", fontsize=12)
    ax.grid(color="0.92", lw=0.6)
    ax.set_axisbelow(True)
    if args.title:
        ax.set_title(args.title, loc="left", fontsize=14)
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, shrink=0.4, anchor=(0, 1))
    cb.set_label("p-adjusted")
    ref = [ratios.min(), ratios.max()]
    handles = [plt.scatter([], [], s=30 + 250 * (v / ratios.max()), color="k", label=f"{v:.2f}") for v in ref]
    ax.legend(handles=handles, title="Gene ratio", loc="upper left", bbox_to_anchor=(1.03, 0.45), frameon=False)
    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(f"{args.out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"Wrote {args.out_prefix}.png/.svg/.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
