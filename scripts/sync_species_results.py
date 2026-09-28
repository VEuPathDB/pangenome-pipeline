#!/usr/bin/env python3
"""
Repo-maintenance tool (runs locally, not on the cluster - nothing here needs SLURM). Builds
results/<Species>/ inside this repository from a species project folder (e.g. ../Afum, ../Tgondii
relative to this repo) so the reports can be published without publishing that whole working
folder, which is full of large/raw intermediate data (hundreds of MB to several GB, including
individual files well over GitHub's 100MB hard limit).

What it does, for one species:
  1. Copies README.md and SUMMARY.md.
  2. Finds every local relative link in those two files (images, tables, other local files) and
     copies each one into results/<Species>/ at the same relative path, so the links keep working
     unchanged. Anything not found under the species folder is reported, not silently skipped.
  3. A link that starts with "../pangenome_analysis/" (this project's own convention for a species
     README linking back to the general pipeline docs) is rewritten to "../../", since
     results/<Species>/ sits two levels below this repo's root instead of one - everything after
     that prefix is left untouched. Links to steps 05-08 will 404 until those steps are published
     here too; that is expected, not a bug in this script.
  4. A "Files" table row that links an entire step folder to itself (this project's SUMMARY.md
     convention, e.g. "| [01_Genome_Selection_QC/](01_Genome_Selection_QC/) | ... |") is de-linked
     to plain text with a note, since that whole folder (raw/intermediate data) is not published -
     only the specific assets actually referenced elsewhere in the text are.
  5. External links (http/https) and in-page anchors (#...) are left exactly as they are.

Nothing outside results/<Species>/ is ever touched, and the species' own working folder
(../Afum, ../Tgondii) is never modified - this only ever reads from it.

Re-run any time the species README/SUMMARY or its referenced assets change; it always rewrites
results/<Species>/ from scratch (removing files that are no longer referenced).

Usage:
    python3 scripts/sync_species_results.py --species Afum
    python3 scripts/sync_species_results.py --species Tgondii --source ../Tgondii --out results/Tgondii
"""

import argparse
import os
import re
import shutil
import sys

LINK_RE = re.compile(r"\]\(([^)]+)\)")
STEP_FOLDER_ROW_RE = re.compile(r"\[(\d\d_[A-Za-z0-9_-]+/)\]\(\1\)")
PIPELINE_LINK_PREFIX = "../pangenome_analysis/"
REPORT_FILES = ("README.md", "SUMMARY.md")


def classify(link):
    """-> one of "external", "anchor", "pipeline", "asset" """
    if link.startswith("http://") or link.startswith("https://") or link.startswith("mailto:"):
        return "external"
    if link.startswith("#"):
        return "anchor"
    if link.startswith(PIPELINE_LINK_PREFIX):
        return "pipeline"
    return "asset"


def rewrite_text(text, source_dir, assets_found, assets_missing):
    # 1. de-link whole-step-folder "Files" table rows - do this before the generic asset pass,
    #    since otherwise those same paths would be treated as (missing, folder-shaped) assets.
    def delink_step_row(m):
        folder = m.group(1)
        return f"{folder} *(working data, not published here)*"

    text = STEP_FOLDER_ROW_RE.sub(delink_step_row, text)

    # 2. every remaining local link: rewrite pipeline cross-repo links, copy local assets in place.
    def handle_link(m):
        link = m.group(1)
        kind = classify(link)
        if kind in ("external", "anchor"):
            return m.group(0)
        if kind == "pipeline":
            return f"]({'../../' + link[len(PIPELINE_LINK_PREFIX) :]})"
        # kind == "asset": a path relative to the species folder, optionally with a "#anchor" tail
        path_part, _, anchor = link.partition("#")
        if path_part in REPORT_FILES or path_part == "":
            return m.group(0)  # README.md/SUMMARY.md linking to each other - already siblings
        src = os.path.join(source_dir, path_part)
        if not os.path.isfile(src):
            assets_missing.append(path_part)
            return m.group(0)
        assets_found.add(path_part)
        return m.group(0)  # path itself is unchanged - only pipeline links change

    text = LINK_RE.sub(handle_link, text)
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True, help="Species name, used as the results/<Species>/ folder name")
    ap.add_argument(
        "--source", help="Path to the species project folder (default: ../<species>, relative to this repo's root)"
    )
    ap.add_argument(
        "--out", help="Where to write the curated copy (default: results/<species>, relative to this repo's root)"
    )
    args = ap.parse_args()

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    source_dir = os.path.abspath(args.source or os.path.join(repo_root, "..", args.species))
    out_dir = os.path.abspath(args.out or os.path.join(repo_root, "results", args.species))

    for f in REPORT_FILES:
        if not os.path.isfile(os.path.join(source_dir, f)):
            sys.exit(f"ERROR: {os.path.join(source_dir, f)} not found")

    if os.path.isdir(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)

    assets_found, assets_missing = set(), []
    texts = {}
    for f in REPORT_FILES:
        with open(os.path.join(source_dir, f), encoding="utf-8") as fh:
            texts[f] = rewrite_text(fh.read(), source_dir, assets_found, assets_missing)

    for f, text in texts.items():
        with open(os.path.join(out_dir, f), "w", encoding="utf-8") as fh:
            fh.write(text)

    copied, copy_failed = [], []
    for rel in sorted(assets_found):
        src = os.path.join(source_dir, rel)
        dst = os.path.join(out_dir, rel)
        try:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(src, dst)
            copied.append(rel)
        except OSError as e:
            copy_failed.append((rel, str(e)))

    total_bytes = sum(os.path.getsize(os.path.join(out_dir, rel)) for rel in copied)
    print(f"[{args.species}] {out_dir}")
    print("  README.md + SUMMARY.md copied")
    print(f"  {len(copied)} asset(s) copied, {total_bytes / 1024 / 1024:.1f} MB total")
    if assets_missing:
        print(
            f"  WARNING: {len(assets_missing)} linked path(s) not found under {source_dir}, left as-is in the text:",
            file=sys.stderr,
        )
        for m in assets_missing:
            print(f"    {m}", file=sys.stderr)
    if copy_failed:
        print(f"  ERROR: {len(copy_failed)} asset(s) failed to copy:", file=sys.stderr)
        for rel, err in copy_failed:
            print(f"    {rel}: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
