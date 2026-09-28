#!/usr/bin/env python3
"""
Summarize TE content across a set of RepeatMasker .out files.

For each genome (one RM_<genome_id>/ or R_<genome_id>/ directory per genome,
as produced by run_repeatmasker.slurm.sh), this script:

  1. Filters the RepeatMasker .out file, dropping Simple_repeat and
     Low_complexity hits, and writes the filtered table to
     <outdir>/filtered_out/<genome_id>.filtered.out
  2. Merges overlapping hits per contig (so a genomic base covered by more
     than one overlapping RepeatMasker hit is only counted once) to compute,
     per genome:
       - total TE bp and % of genome
       - bp and % per superfamily (RepeatMasker "class/family" column)
       - bp and % per family (RepeatMasker "matching repeat" column)
  3. Writes summary TSVs and a set of publication-quality plots (matplotlib, static PNG).

Usage:
    python3 -m pangenome_pipeline.step02_te_annotation_clean_up.summarize_te_content --indir . \
        --outdir te_content_summary

Expects, per genome, a directory named RM_<genome_id>/ (or R_<genome_id>/) under --indir
containing exactly one RepeatMasker .out file and one .tbl file. Genome length is read from
the .tbl file's "total length" line.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Validated categorical palette (dataviz skill reference palette, light mode).
# Fixed order - never re-cycled/re-sorted per chart.
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
# Extended palette for the family-composition plot, which can have more
# individually-shown categories than the 8-hue base palette covers. The
# first 8 entries are identical to CATEGORICAL (same hues, same order) so a
# family that also happens to be a top superfamily still reads consistently;
# the next 8 extend it so up to 16 families can get a distinct color before
# folding into "Other".
EXTENDED_CATEGORICAL = CATEGORICAL + [
    "#8ecae6",
    "#bc6c25",
    "#6a4c93",
    "#606c38",
    "#ffb703",
    "#457b9d",
    "#d62828",
    "#a8dadc",
]
OTHER_COLOR = "#c3c2b7"  # neutral, recessive - for the folded "Other" bucket

SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"

REMOVE_CLASSES = {"simple_repeat", "low_complexity"}

TBL_LENGTH_RE = re.compile(r"total length:\s*([\d,]+)\s*bp", re.IGNORECASE)


def find_genome_dirs(indir: Path):
    """RM_<genome_id> or R_<genome_id> subdirectories, each with one .out and one .tbl file."""
    genomes = {}
    for d in sorted(indir.glob("R*_*")):
        if not d.is_dir():
            continue
        if d.name.startswith("RM_"):
            genome_id = d.name[len("RM_") :]
        elif d.name.startswith("R_"):
            genome_id = d.name[len("R_") :]
        else:
            genome_id = d.name
        out_files = list(d.glob("*.out"))
        tbl_files = list(d.glob("*.tbl"))
        if not out_files or not tbl_files:
            print(f"WARNING: skipping {d} (missing .out or .tbl)")
            continue
        genomes[genome_id] = {"out": out_files[0], "tbl": tbl_files[0]}
    return genomes


def genome_length_from_tbl(tbl_path: Path) -> int:
    text = tbl_path.read_text()
    m = TBL_LENGTH_RE.search(text)
    if not m:
        raise ValueError(f"Could not find 'total length' in {tbl_path}")
    return int(m.group(1).replace(",", ""))


def parse_out(out_path: Path):
    """Yield (seq, start, end, repeat_name, repeat_class, seq_left, pct_div,
    rm_id) for every hit, including Simple_repeat/Low_complexity (filtering
    happens in the caller so we can still write out the pre-filter line
    count if ever needed). seq_left is RepeatMasker's "(left)" column for
    the query sequence - the number of bases in the contig past this hit's
    end - used to infer each contig's total length
    (contig_length = end + seq_left). pct_div is RepeatMasker's own "%
    divergence from consensus" column. rm_id is RepeatMasker's own ID
    column (last field): fragments of one element interrupted by a nested
    insertion of another element share the same ID, on the same contig."""
    hits = []
    with out_path.open() as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            fields = line.split()
            if len(fields) < 15:
                continue
            seq = fields[4]
            start = int(fields[5])
            end = int(fields[6])
            seq_left = int(fields[7].strip("()"))
            pct_div = float(fields[1])
            # Columns are fixed regardless of strand (+ or C):
            # repeat name is always field 10, class/family is always field 11.
            repeat_name = fields[9]
            repeat_class = fields[10]
            rm_id = fields[14]
            hits.append((seq, start, end, repeat_name, repeat_class, seq_left, pct_div, rm_id))
    return hits


def contig_lengths(hits) -> dict:
    """Infer each contig's total length from RepeatMasker's own "(left)"
    column (contig_length = hit_end + bases_left_in_query), using every hit
    on that contig (any hit gives the same answer; using them all is just
    robustness against a malformed line). Contigs with zero RepeatMasker
    hits at all are invisible to the .out file and cannot be recovered here."""
    lengths: dict[str, int] = {}
    for seq, start, end, name, cls, seq_left, pct_div, rm_id in hits:
        length = end + seq_left
        if length > lengths.get(seq, 0):
            lengths[seq] = length
    return lengths


def merge_bp(intervals_by_seq):
    """intervals_by_seq: dict seq -> list of (start,end), 1-based inclusive.
    Returns total non-redundant bp across all sequences."""
    total = 0
    for intervals in intervals_by_seq.values():
        if not intervals:
            continue
        intervals = sorted(intervals)
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e + 1:
                cur_e = max(cur_e, e)
            else:
                total += cur_e - cur_s + 1
                cur_s, cur_e = s, e
        total += cur_e - cur_s + 1
    return total


def write_filtered_out(out_path: Path, hits, dest: Path):
    """Re-write the original .out file, dropping filtered-class lines, and
    keeping the original 3-line header intact for readability."""
    header_lines = []
    with out_path.open() as fh:
        for i, line in enumerate(fh):
            if i < 3:
                header_lines.append(line)
    kept_lines = []
    with out_path.open() as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            fields = line.split()
            if len(fields) < 14:
                continue
            repeat_class = fields[10]
            if repeat_class.lower() in REMOVE_CLASSES:
                continue
            kept_lines.append(line)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w") as out:
        out.writelines(header_lines)
        out.writelines(kept_lines)


def fold_other(series: pd.Series, top_n: int, palette=CATEGORICAL) -> pd.Series:
    """Keep the top_n categories by value, fold the rest into 'Other'.
    top_n is capped at the size of palette so a chart never has to repeat a
    hue for two different categories."""
    top_n = min(top_n, len(palette))
    if len(series) <= top_n:
        return series.sort_values(ascending=False)
    top = series.sort_values(ascending=False).iloc[:top_n]
    other = series.sort_values(ascending=False).iloc[top_n:].sum()
    result = top.copy()
    result["Other"] = other
    return result


def style_axes(ax, horizontal=False):
    ax.set_facecolor(SURFACE)
    ax.figure.set_facecolor(SURFACE)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(BASELINE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    if horizontal:
        ax.xaxis.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    else:
        ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.title.set_color(INK_PRIMARY)
    ax.xaxis.label.set_color(INK_SECONDARY)
    ax.yaxis.label.set_color(INK_SECONDARY)


def plot_bar(values: pd.Series, title: str, xlabel: str, dest: Path, color=CATEGORICAL[0]):
    """Horizontal bar chart, one bar per genome, largest at the top. Figure
    height tracks the number of genomes/strains directly (a small constant
    plus a fixed increment per bar) so the gap between bars stays the same
    regardless of how many strains are plotted, instead of growing with a
    generous floor meant for small genome sets."""
    values = values.sort_values()
    fig, ax = plt.subplots(figsize=(7, 0.28 * len(values) + 1.2))
    ax.barh(values.index, values.to_numpy(), color=color, height=0.65, zorder=3)
    style_axes(ax, horizontal=True)
    ax.ticklabel_format(style="plain", axis="x")
    ax.xaxis.set_major_formatter(lambda x, pos: f"{x:,.0f}")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_scatter(
    df: pd.DataFrame, xcol: str, ycol: str, title: str, xlabel: str, ylabel: str, dest: Path, color=CATEGORICAL[0]
):
    """One point per genome, labeled - for spotting a relationship between
    two per-genome metrics (e.g. does more fragmentation/divergence line up
    with less total TE detected?)."""
    df = df.dropna(subset=[xcol, ycol])
    if df.empty:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.scatter(df[xcol], df[ycol], color=color, s=50, zorder=3, edgecolor=SURFACE, linewidth=0.8)
    # stagger the label offsets a little so points in a tight cluster don't
    # write over each other
    for i, (genome_id, row) in enumerate(df.sort_values(xcol).iterrows()):
        dy = 6 if i % 2 == 0 else -10
        ax.annotate(
            str(genome_id),
            (row[xcol], row[ycol]),
            fontsize=7,
            color=INK_SECONDARY,
            xytext=(5, dy),
            textcoords="offset points",
        )
    style_axes(ax, horizontal=True)
    ax.yaxis.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    # let matplotlib pick the tick labels: these axes are narrow-range decimals
    # (% divergence, % genome), not the large bp counts the bar charts format
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_stacked_pct(df: pd.DataFrame, title: str, dest: Path, palette=CATEGORICAL):
    """df: rows = genome, columns = category (already folded to <=len(palette)
    + Other), values = % of genome. Horizontal 100%-style stacked bar, one
    row per genome, largest total at the top."""
    categories = list(df.columns)
    colors = {}
    cat_no_other = [c for c in categories if c != "Other"]
    for i, cat in enumerate(cat_no_other):
        colors[cat] = palette[i % len(palette)]
    if "Other" in categories:
        colors["Other"] = OTHER_COLOR

    df = df.loc[df.sum(axis=1).sort_values().index]

    fig, ax = plt.subplots(figsize=(8, max(4, 0.45 * len(df) + 1.5)))
    left = pd.Series([0.0] * len(df), index=df.index)
    for cat in categories:
        ax.barh(
            df.index,
            df[cat],
            left=left,
            label=cat,
            color=colors[cat],
            height=0.6,
            zorder=3,
            edgecolor=SURFACE,
            linewidth=1,
        )
        left = left + df[cat]
    style_axes(ax, horizontal=True)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.set_xlabel("% of genome")
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8, labelcolor=INK_SECONDARY)
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--indir",
        required=True,
        help="Directory containing RM_<genome_id>/ or R_<genome_id>/ subdirectories "
        "(each with one .out and one .tbl file)",
    )
    ap.add_argument("--outdir", required=True, help="Directory to write filtered .out files, summary TSVs, and plots")
    ap.add_argument(
        "--top-superfamilies",
        type=int,
        default=8,
        help="Number of superfamilies to show individually in "
        "the composition plot before folding the rest into "
        "'Other' (default: 8)",
    )
    ap.add_argument(
        "--top-families",
        type=int,
        default=15,
        help="Number of families to show individually in the "
        "family composition plot before folding the rest "
        "into 'Other' (default: 15, capped at the size of "
        "the extended color palette so hues never repeat - "
        "the full, unfolded per-family data is always in "
        "family_summary.tsv)",
    )
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    filtered_dir = outdir / "filtered_out"
    plots_dir = outdir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    genomes = find_genome_dirs(indir)
    if not genomes:
        raise SystemExit(f"ERROR: no RM_<genome_id>/ or R_<genome_id>/ directories with .out+.tbl found in {indir}")

    genome_rows = []
    superfamily_rows = []
    family_rows = []

    for genome_id, paths in sorted(genomes.items()):
        print(f"Processing {genome_id} ...")
        genome_len = genome_length_from_tbl(paths["tbl"])
        all_hits = parse_out(paths["out"])

        write_filtered_out(paths["out"], all_hits, filtered_dir / f"{genome_id}.filtered.out")

        kept_hits = [h for h in all_hits if h[4].lower() not in REMOVE_CLASSES]

        overall_intervals = defaultdict(list)
        by_class_intervals = defaultdict(lambda: defaultdict(list))
        by_family_intervals = defaultdict(lambda: defaultdict(list))
        div_weighted_sum = 0.0
        div_weight_total = 0.0
        family_div_weighted = defaultdict(float)
        family_div_weight = defaultdict(float)
        # Group hit lines by (contig, RepeatMasker ID): pieces of one element
        # interrupted by a nested insertion of another element share the same
        # ID on the same contig. Treating each piece as its own "fragment"
        # would inflate apparent fragmentation from nested insertions alone -
        # a confound unrelated to divergence from the curated library, and
        # roughly constant across genomes. Grouping first, then taking each
        # group's full span (first piece's start to last piece's end),
        # recovers the true size of the original inserted copy.
        element_spans = {}
        for seq, start, end, name, cls, seq_left, pct_div, rm_id in kept_hits:
            overall_intervals[seq].append((start, end))
            by_class_intervals[cls][seq].append((start, end))
            by_family_intervals[name][seq].append((start, end))
            hit_len = end - start + 1
            div_weighted_sum += pct_div * hit_len
            div_weight_total += hit_len
            family_div_weighted[name] += pct_div * hit_len
            family_div_weight[name] += hit_len

            key = (seq, rm_id)
            lo, hi, _ = element_spans.get(key, (start, end, cls))
            element_spans[key] = (min(lo, start), max(hi, end), cls)

        num_fragments = len(element_spans)
        mean_pct_divergence = div_weighted_sum / div_weight_total if div_weight_total else None

        contig_lens = contig_lengths(all_hits)
        contig_rows = []
        for seq, intervals in overall_intervals.items():
            contig_bp = merge_bp({seq: intervals})
            contig_len = contig_lens.get(seq)
            contig_rows.append(
                {
                    "genome_id": genome_id,
                    "contig": seq,
                    "contig_length_bp": contig_len,
                    "TE_bp": contig_bp,
                    "TE_pct": 100.0 * contig_bp / contig_len if contig_len else None,
                }
            )
        contig_df = pd.DataFrame(contig_rows).sort_values("contig_length_bp", ascending=False)
        contig_dir = outdir / "per_contig"
        contig_dir.mkdir(parents=True, exist_ok=True)
        contig_df.to_csv(contig_dir / f"{genome_id}.te_per_contig.tsv", sep="\t", index=False)

        total_te_bp = merge_bp(overall_intervals)
        total_te_pct = 100.0 * total_te_bp / genome_len if genome_len else 0.0
        genome_rows.append(
            {
                "genome_id": genome_id,
                "genome_length_bp": genome_len,
                "total_TE_bp": total_te_bp,
                "total_TE_pct": total_te_pct,
                "num_TE_fragments": num_fragments,
                "avg_fragment_length_bp": total_te_bp / num_fragments if num_fragments else None,
                "mean_pct_divergence": mean_pct_divergence,
            }
        )

        for cls, intervals_by_seq in by_class_intervals.items():
            bp = merge_bp(intervals_by_seq)
            superfamily_rows.append(
                {
                    "genome_id": genome_id,
                    "superfamily": cls,
                    "bp": bp,
                    "pct_of_genome": 100.0 * bp / genome_len if genome_len else 0.0,
                }
            )

        for name, intervals_by_seq in by_family_intervals.items():
            bp = merge_bp(intervals_by_seq)
            weight = family_div_weight[name]
            family_rows.append(
                {
                    "genome_id": genome_id,
                    "family": name,
                    "bp": bp,
                    "pct_of_genome": 100.0 * bp / genome_len if genome_len else 0.0,
                    "mean_pct_divergence": family_div_weighted[name] / weight if weight else None,
                }
            )

    genome_df = pd.DataFrame(genome_rows).set_index("genome_id").sort_index()
    superfamily_df = pd.DataFrame(superfamily_rows)
    family_df = pd.DataFrame(family_rows)

    genome_df.to_csv(outdir / "te_content_summary.tsv", sep="\t")
    superfamily_df.to_csv(outdir / "superfamily_summary.tsv", sep="\t", index=False)
    family_df.to_csv(outdir / "family_summary.tsv", sep="\t", index=False)
    print(f"\nWrote {outdir / 'te_content_summary.tsv'}")
    print(f"Wrote {outdir / 'superfamily_summary.tsv'}")
    print(f"Wrote {outdir / 'family_summary.tsv'}")

    # --- Plots ---------------------------------------------------------
    # bp values are shown in kb on the plots (avoids scientific notation on
    # the axis for whole-genome sizes); the underlying TSVs stay in bp.
    plot_bar(
        genome_df["genome_length_bp"] / 1000.0,
        "Genome size (kb)",
        "kb",
        plots_dir / "genome_size_bp.png",
        color=CATEGORICAL[2],
    )
    plot_bar(genome_df["total_TE_bp"] / 1000.0, "Total TE content (kb)", "kb", plots_dir / "total_te_bp.png")
    plot_bar(
        genome_df["total_TE_pct"],
        "Total TE content (% of genome)",
        "%",
        plots_dir / "total_te_pct.png",
        color=CATEGORICAL[1],
    )

    # Tests whether divergence from the (single-strain-based) curated
    # library is under-detecting TEs in more divergent genomes: higher mean
    # divergence lining up with less total TE content would support that.
    plot_scatter(
        genome_df,
        "mean_pct_divergence",
        "total_TE_pct",
        "Mean divergence from consensus vs. total TE content",
        "Mean % divergence from consensus",
        "Total TE (% of genome)",
        plots_dir / "divergence_vs_total_te.png",
        color=CATEGORICAL[4],
    )

    superfamily_pivot = superfamily_df.pivot_table(
        index="genome_id", columns="superfamily", values="pct_of_genome", fill_value=0.0
    )
    superfamily_totals = superfamily_pivot.sum(axis=0)
    folded_categories = fold_other(superfamily_totals, args.top_superfamilies).index.tolist()
    plot_cols = [c for c in folded_categories if c != "Other"]
    superfamily_plot_df = pd.DataFrame(index=superfamily_pivot.index)
    for c in plot_cols:
        superfamily_plot_df[c] = superfamily_pivot.get(c, 0.0)
    other_cols = [c for c in superfamily_pivot.columns if c not in plot_cols]
    if other_cols:
        superfamily_plot_df["Other"] = superfamily_pivot[other_cols].sum(axis=1)
    superfamily_plot_df.to_csv(outdir / "superfamily_composition_table.tsv", sep="\t")
    plot_stacked_pct(
        superfamily_plot_df, "TE superfamily composition (% of genome)", plots_dir / "superfamily_composition.png"
    )

    family_pivot = family_df.pivot_table(index="genome_id", columns="family", values="pct_of_genome", fill_value=0.0)
    family_totals = family_pivot.sum(axis=0)
    folded_fam_categories = fold_other(family_totals, args.top_families, palette=EXTENDED_CATEGORICAL).index.tolist()
    fam_plot_cols = [c for c in folded_fam_categories if c != "Other"]
    family_plot_df = pd.DataFrame(index=family_pivot.index)
    for c in fam_plot_cols:
        family_plot_df[c] = family_pivot.get(c, 0.0)
    other_fam_cols = [c for c in family_pivot.columns if c not in fam_plot_cols]
    if other_fam_cols:
        family_plot_df["Other"] = family_pivot[other_fam_cols].sum(axis=1)
    family_plot_df.to_csv(outdir / "family_composition_table.tsv", sep="\t")
    plot_stacked_pct(
        family_plot_df,
        "TE family composition (% of genome, top families)",
        plots_dir / "family_composition.png",
        palette=EXTENDED_CATEGORICAL,
    )

    print(f"\nPlots written to {plots_dir}/")
    print(
        f"Composition tables (matching the plots) written to {outdir}/"
        "superfamily_composition_table.tsv and family_composition_table.tsv"
    )


if __name__ == "__main__":
    main()
