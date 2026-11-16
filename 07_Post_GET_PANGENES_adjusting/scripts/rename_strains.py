#!/usr/bin/env python3
"""
Show the strains under other names (e.g. GCA_040126065.1 -> AF100-9B) in the files made for people: the explorer
and the plots. The pipeline itself keeps the original names (the step 01-07 folders, GFFs and tables use them).

The names table has two columns, old name and new name, tab-separated; lines starting with # and a header line
"old<TAB>new" are skipped. New names must be unique and should have no spaces.

Two modes:
  --tsv IN OUT       a matrix or table: every cell equal to an old name, or starting with "<old name>_" (get_pangenes
                     headers such as <genome>_for_get_pangenes...), gets the new name.
  --json-dir DIR     the explorer's data folder, changed in place: every JSON string equal to an old name (as a value
                     or a key) gets the new name, and every file <prefix>_<old name>.json is renamed <prefix>_<new name>.json.
  --name NAME        print the new name of NAME (NAME itself when it is not in the table); several names separated
                     by commas are mapped one by one.

Usage:
    python3 rename_strains.py --names strain_names.tsv --tsv pangene_matrix_genes.tr.tab matrix_named.tr.tab
    python3 rename_strains.py --names strain_names.tsv --json-dir pangenome_explorer/final/data
    python3 rename_strains.py --names strain_names.tsv --name AF293

Needs only the standard library.
"""
import argparse
import csv
import glob
import json
import os
import sys


def read_names(path):
    names = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = [x.strip() for x in line.rstrip("\n").split("\t")]
            if len(f) < 2 or not f[0] or not f[1]:
                sys.exit(f"ERROR: {path}: needs two tab-separated columns (old, new): {line!r}")
            if f[0].lower() == "old" and f[1].lower() == "new":
                continue
            if f[0] in names:
                sys.exit(f"ERROR: {path}: {f[0]} is listed twice")
            names[f[0]] = f[1]
    news = list(names.values())
    dup = sorted({n for n in news if news.count(n) > 1})
    if dup:
        sys.exit(f"ERROR: {path}: new names used twice: {', '.join(dup)}")
    clash = sorted(n for n in news if n in names and names[n] != n)
    if clash:
        sys.exit(f"ERROR: {path}: new names that are also old names of other strains: {', '.join(clash)}")
    return {o: n for o, n in names.items() if o != n}


def rename_cell(c, names):
    if c in names:
        return names[c]
    for o, n in names.items():
        if c.startswith(o + "_"):
            return n + c[len(o):]
    return c


def do_tsv(src, out, names):
    with open(src, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    n = 0
    for r in rows:
        for i, c in enumerate(r):
            new = rename_cell(c, names)
            if new != c:
                r[i] = new
                n += 1
    with open(out, "w", newline="") as fo:
        csv.writer(fo, delimiter="\t", lineterminator="\n").writerows(rows)
    print(f"{src} -> {out}: {n} cells renamed", file=sys.stderr)


def walk(x, names, cnt):
    if isinstance(x, str):
        if x in names:
            cnt[0] += 1
            return names[x]
        return x
    if isinstance(x, list):
        return [walk(v, names, cnt) for v in x]
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            if k in names:
                cnt[0] += 1
                k = names[k]
            out[k] = walk(v, names, cnt)
        return out
    return x


def do_json_dir(d, names):
    files = sorted(glob.glob(os.path.join(d, "*.json")))
    if not files:
        sys.exit(f"ERROR: no JSON files in {d}")
    tot = 0
    for p in files:
        with open(p) as fh:
            data = json.load(fh)
        cnt = [0]
        data = walk(data, names, cnt)
        base = os.path.basename(p)[:-5]
        new_base = base
        for o, n in names.items():
            if base.endswith("_" + o):
                new_base = base[: -len(o)] + n
                break
        new_p = os.path.join(d, new_base + ".json")
        # a file with the new name left by an earlier run is replaced (the builders rewrite the old names each run)
        with open(new_p, "w") as fo:
            json.dump(data, fo, separators=(",", ":"))
        if new_p != p:
            os.remove(p)
        tot += cnt[0]
        if cnt[0] or new_p != p:
            print(f"{base}.json" + (f" -> {new_base}.json" if new_p != p else "") + f": {cnt[0]} names", file=sys.stderr)
    print(f"{d}: {len(files)} files, {tot} names changed", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--names", required=True, help="table: old name <TAB> new name")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--tsv", nargs=2, metavar=("IN", "OUT"))
    g.add_argument("--json-dir")
    g.add_argument("--name")
    args = ap.parse_args()
    names = read_names(args.names)
    if args.tsv:
        do_tsv(args.tsv[0], args.tsv[1], names)
    elif args.json_dir:
        do_json_dir(args.json_dir, names)
    else:
        print(",".join(names.get(x, x) for x in args.name.split(",")) if args.name else "")


if __name__ == "__main__":
    main()
