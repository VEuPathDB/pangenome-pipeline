#!/usr/bin/env python3
"""
Fix FASTQ read headers for Trinity / funannotate compatibility.

Requirements for Trinity compatibility:
  - Header line: '@<accession>.<read_id> length=<len>' (or simple ID)
  - 3rd line: Must be exactly '+' without the header repeat.
"""
import gzip
import os
import sys

# =============================================================================
# CONFIGURATION — edit these paths for your run
# =============================================================================
LEFT_INPUT   = 'combined_left_trinity.fastq.gz'
RIGHT_INPUT  = 'combined_right_trinity.fastq.gz'

LEFT_OUTPUT  = 'combined_left_trinity.fixed_headers.fastq.gz'
RIGHT_OUTPUT = 'combined_right_trinity.fixed_headers.fastq.gz'
# =============================================================================


def fix_header(header_line, suffix):
    header_line = header_line.rstrip('\r\n')
    if header_line.startswith('@'):
        header_line = header_line[1:]
    parts = header_line.split()
    read_id = parts[0]
    return f"@{read_id}{suffix}"


def process_file(in_path, out_path, suffix):
    print(f"Processing: {in_path} -> {out_path}")
    count = 0
    with gzip.open(in_path, 'rt') as fin, gzip.open(out_path, 'wt') as fout:
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            fixed_header = fix_header(header, suffix)
            fout.write(fixed_header + '\n')
            fout.write(seq)
            fout.write('+\n')
            fout.write(qual)

            count += 1
            if count % 10000000 == 0:
                print(f"  Processed {count:,} reads...")
    print(f"Finished {out_path}: {count:,} reads total.")


def main():
    left_in = sys.argv[1] if len(sys.argv) > 1 else LEFT_INPUT
    right_in = sys.argv[2] if len(sys.argv) > 2 else RIGHT_INPUT
    left_out = sys.argv[3] if len(sys.argv) > 3 else LEFT_OUTPUT
    right_out = sys.argv[4] if len(sys.argv) > 4 else RIGHT_OUTPUT

    if not os.path.exists(left_in):
        print(f"Error: {left_in} not found.")
        sys.exit(1)
    if not os.path.exists(right_in):
        print(f"Error: {right_in} not found.")
        sys.exit(1)

    process_file(left_in, left_out, "/1")
    process_file(right_in, right_out, "/2")
    print("All headers fixed successfully.")


if __name__ == '__main__':
    main()
