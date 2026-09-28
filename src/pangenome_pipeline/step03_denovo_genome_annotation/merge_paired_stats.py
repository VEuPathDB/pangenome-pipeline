#!/usr/bin/env python3
"""
merge_paired_stats.py

Merge Trinity left/right k-mer statistics when the default Trinity Perl merge
fails to match paired SRA read IDs.

Expected input files in the current working directory:
    left.fa.K25.stats.sort
    right.fa.K25.stats.sort

Output:
    pairs.K25.stats.custom

This script removes the mate suffixes "_1/1" and "_2/2" before matching read
IDs, then writes the paired stats file expected by Trinity's nbkc_normalize.pl.
"""

import sys
from pathlib import Path

LEFT_STATS = Path("left.fa.K25.stats.sort")
RIGHT_STATS = Path("right.fa.K25.stats.sort")
OUTPUT_STATS = Path("pairs.K25.stats.custom")


def check_input(path):
    """Fail early if a required input file is missing or empty."""
    if not path.exists():
        sys.exit(f"ERROR: Missing required file: {path}")
    if path.stat().st_size == 0:
        sys.exit(f"ERROR: Required file is empty: {path}")


def load_stats(path, suffix):
    """
    Load one Trinity stats file into a dictionary keyed by the base read ID.
    Example:
        SRR32610147.1_1/1 -> SRR32610147.1
        SRR32610147.1_2/2 -> SRR32610147.1
    """
    stats = {}
    malformed = 0

    with path.open() as handle:
        header = next(handle, None)
        if header is None:
            sys.exit(f"ERROR: Empty stats file: {path}")

        for line_number, line in enumerate(handle, start=2):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                malformed += 1
                continue

            base_id = parts[0].replace(suffix, "")
            stats[base_id] = parts

            if len(stats) % 5_000_000 == 0:
                print(f"Loaded {len(stats):,} records from {path}", flush=True)

    if malformed:
        print(f"WARNING: Skipped {malformed:,} malformed lines in {path}", file=sys.stderr, flush=True)

    print(f"Loaded {len(stats):,} usable records from {path}", flush=True)
    return stats


def merge_paired_stats(left_file, right_file, output_file):
    """Merge left and right Trinity stats into one paired stats file."""
    check_input(left_file)
    check_input(right_file)

    print("Loading left stats...", flush=True)
    left_stats = load_stats(left_file, "_1/1")

    print("Loading right stats...", flush=True)
    right_stats = load_stats(right_file, "_2/2")

    print("Merging paired stats...", flush=True)
    pairs_written = 0
    left_without_right = 0

    with output_file.open("w") as output_handle:
        output_handle.write(
            "acc\tleft_acc\tleft_median_cov\tleft_mean_cov\tleft_stdev\t"
            "right_acc\tright_median_cov\tright_mean_cov\tright_stdev\t"
            "median_cov\tmean_cov\tstdev\n"
        )

        for base_id in sorted(left_stats):
            right = right_stats.get(base_id)
            if right is None:
                left_without_right += 1
                continue

            left = left_stats[base_id]

            try:
                median_cov = (float(left[1]) + float(right[1])) / 2
                mean_cov = (float(left[2]) + float(right[2])) / 2
                stdev = (float(left[3]) + float(right[3])) / 2
            except (IndexError, ValueError) as error:
                print(f"WARNING: Skipping malformed pair {base_id}: {error}", file=sys.stderr, flush=True)
                continue

            output_handle.write(
                f"{base_id}\t"
                f"{left[0]}\t{left[1]}\t{left[2]}\t{left[3]}\t"
                f"{right[0]}\t{right[1]}\t{right[2]}\t{right[3]}\t"
                f"{median_cov}\t{mean_cov}\t{stdev}\n"
            )

            pairs_written += 1
            if pairs_written % 5_000_000 == 0:
                print(f"Merged {pairs_written:,} paired reads", flush=True)

    if pairs_written == 0:
        sys.exit(
            "ERROR: No paired reads were merged. Check whether the read IDs use the expected suffixes _1/1 and _2/2."
        )

    print(f"Successfully merged {pairs_written:,} paired reads", flush=True)
    if left_without_right:
        print(f"WARNING: {left_without_right:,} left reads had no matching right read", file=sys.stderr, flush=True)
    print(f"Output written to: {output_file}", flush=True)


def main():
    merge_paired_stats(LEFT_STATS, RIGHT_STATS, OUTPUT_STATS)


if __name__ == "__main__":
    main()
