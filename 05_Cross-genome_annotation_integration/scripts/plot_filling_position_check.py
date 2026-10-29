"""Scheme: filling a genome with a Liftoff transfer, and the check of its neighbours. Illustrative genes.
Usage: python3 plot_filling_position_check.py <out.png>"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

COL = {"A": "#9AA5AE", "B": "#9AA5AE", "C": "#9AA5AE", "D": "#9AA5AE", "P": "#C9B27C", "Q": "#C9B27C",
       "R": "#C9B27C", "S": "#C9B27C", "X": "#4C72B0", "X'": "#8172B2"}


def chrom(ax, y, x0, genes, label, title=None, empty_at=None):
    ax.plot([x0 - 0.3, x0 + len(genes) * 1.3 + 0.1], [y, y], color="#C9D1D6", lw=2, zorder=0)
    ax.text(x0 - 0.5, y, label, ha="right", va="center", fontsize=11)
    for i, g in enumerate(genes):
        x = x0 + i * 1.3
        if g is None:
            ax.add_patch(FancyBboxPatch((x, y - 0.22), 1.0, 0.44, boxstyle="round,pad=0.02", fill=False,
                                        ec="#D62728", ls="--", lw=1.5))
            ax.text(x + 0.5, y - 0.45, "nothing\nannotated", ha="center", va="top", fontsize=8, color="#D62728")
            continue
        ax.add_patch(FancyBboxPatch((x, y - 0.22), 1.0, 0.44, boxstyle="round,pad=0.02", fc=COL[g], ec="none"))
        ax.text(x + 0.5, y, g, ha="center", va="center", fontsize=11, color="white", fontweight="bold")
    if title:
        ax.text(x0 - 0.3, y + 0.55, title, fontsize=10, color="#52606D")


fig, ax = plt.subplots(figsize=(13, 8.5))
ax.text(0, 9.6, "1. In most strains, gene X is on chromosome 2, between genes A, B and C, D", fontsize=13,
        fontweight="bold")
chrom(ax, 8.9, 1.5, ["A", "B", "X", "C", "D"], "ARI", "chromosome 2")
chrom(ax, 8.1, 1.5, ["A", "B", "X", "C", "D"], "CAST")
chrom(ax, 7.3, 1.5, ["A", "B", "X", "C", "D"], "…")
ax.text(10.2, 8.55, "X has a similar copy X' elsewhere\n(a paralogue), on chromosome 5:", fontsize=10, va="center")
chrom(ax, 7.3, 10.5, ["P", "Q", "X'", "R", "S"], "", "chromosome 5")

ax.text(0, 6.2, "2. VEG's own annotation has no X. The pipeline takes the Liftoff copy of X and puts it in VEG",
        fontsize=13, fontweight="bold")
chrom(ax, 5.3, 1.5, ["A", "B", None, "C", "D"], "VEG", "chromosome 2")
chrom(ax, 5.3, 10.5, ["P", "Q", "X'", "R", "S"], "", "chromosome 5")
ax.text(4.4, 4.25, "Liftoff places X where its sequence aligns best: usually here (right place) …",
        fontsize=10, color="#2E7D6B", ha="center")
ax.text(13.2, 4.25, "… but sometimes here, over the copy X' (wrong place)", fontsize=10, color="#D62728",
        ha="center")

ax.text(0, 3.3, "3. The check: look at the genes around the place where the copy landed", fontsize=13,
        fontweight="bold")
chrom(ax, 2.4, 1.5, ["A", "B", "X", "C", "D"], "VEG")
ax.text(4.4, 1.6, "neighbours A, B, C, D = the same as in the other strains\n→ right place: X is added",
        ha="center", va="top", fontsize=10, color="#2E7D6B", fontweight="bold")
chrom(ax, 2.4, 10.5, ["P", "Q", "X", "R", "S"], "")
ax.text(13.2, 1.6, "neighbours P, Q, R, S = not the neighbours of X\n→ wrong place: X is NOT added here",
        ha="center", va="top", fontsize=10, color="#D62728", fontweight="bold")
ax.set_xlim(-1.2, 17.5)
ax.set_ylim(0.4, 10)
ax.axis("off")
fig.savefig(sys.argv[1], dpi=170, bbox_inches="tight")
print(sys.argv[1])
