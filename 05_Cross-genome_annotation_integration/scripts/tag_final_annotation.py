#!/usr/bin/env python3
"""
Alternative to filter_final_annotation.py: instead of dropping any gene, tag every gene in
every strain's <genome>_final_annotation.gff with what its CDS overlaps, and pass everything
through unfiltered to the get_pangenes prep step. The idea: some overlaps are local artifacts
of one genome's Liftoff transfer (see the general README's overlap-rule notes) - deleting a
gene before get_pangenes has run risks a false gene-absence call for a real, panel-wide gene
that only misbehaves in one genome. Tagging instead defers that judgment to after get_pangenes,
where the pangene's real presence/absence pattern across the whole panel is known and is much
stronger evidence than anything computable from a single genome's GFF alone.

Adds one new attribute to every protein-coding gene line (not its children):
    intersect=NA                        no CDS-to-CDS overlap >= --min-overlap-bp with anything
    intersect=with_main                 overlaps the reference-derived "main" gene
    intersect=with_other_with_domain    overlaps a non-main gene that itself has a protein domain
    intersect=with_other_no_domain      overlaps a non-main gene with no protein domain either

If a gene overlaps more than one other gene (rare, but possible), the tag is the single most
informative of the above - priority with_main > with_other_with_domain > with_other_no_domain -
not a list of every partner (that detail is still in --out-dir/<genome>_intersect_pairs.tsv,
one row per overlapping pair, for anyone who wants to look closer).

Protein domain = the gene has an InterPro or Pfam hit in its own Dbxref attribute - same
detection as filter_final_annotation.py, see that script's docstring for the caveat on
genomes whose original annotation never recorded this field.

This script never removes a gene. It does, however, do the same two get_pangenes-readiness
things filter_final_annotation.py does, so either script's output is ready to hand straight to
get_pangenes without an extra prep pass:

- every protein-coding gene is relabeled to feature type "gene" (get_pangenes only clusters
  features literally typed that), and its "pseudogenic_transcript" child (if any) to
  "transcript" - non-coding RNA genes are left exactly as they are, feature type and source
  both, since they were already retyped away from "gene" back in step 3
  (prepare_for_liftoff.py) specifically so get_pangenes would never see them as candidates;
- every gene's ID is renumbered to "<id-prefix><strain-index>g<00001, 00002, ...>", same
  scheme and same reason as filter_final_annotation.py (the same source gene ID can otherwise
  appear unchanged in several genomes' GFFs, since this pipeline's Liftoff transfers preserve
  it - get_pangenes needs every ID unique across the whole panel). The original ID is kept in
  a new "original_ID=" attribute.

Usage:
    python3 tag_final_annotation.py \\
        --strains genome_order.txt \\
        --final-dir expanded_annotation \\
        --out-dir tagged_annotation
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from pangenome_pipeline.common.io import read_strains

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotation_common import (
    CHILD_FEATURE_TYPES,
    ID_RE,
    add_original_id_attr,
    assign_new_ids,
    find_cds_overlaps,
    parse_gff,
    replace_attr_value,
)

TAG_PRIORITY = {"with_main": 3, "with_other_with_domain": 2, "with_other_no_domain": 1}


def compute_tags(genes, overlap_pairs, ref_tag):
    """gene_id -> (tag, [(partner_id, overlap_bp), ...]) for every gene involved in >=1 pair."""
    best_tag = {}
    partners = defaultdict(list)

    def update(gid, tag, partner_id, ov_bp):
        partners[gid].append((partner_id, ov_bp))
        cur = best_tag.get(gid)
        if cur is None or TAG_PRIORITY[tag] > TAG_PRIORITY[cur]:
            best_tag[gid] = tag

    for g1_id, g2_id, ov_bp in overlap_pairs:
        g1, g2 = genes[g1_id], genes[g2_id]
        for gid, own, partner in [(g1_id, g1, g2), (g2_id, g2, g1)]:
            if partner["source"] == ref_tag:
                tag = "with_main"
            elif partner["has_domain"]:
                tag = "with_other_with_domain"
            else:
                tag = "with_other_no_domain"
            update(gid, tag, partner["id"], ov_bp)

    return best_tag, partners


def write_tagged_gff(gff_path, out_path, best_tag, id_mapping):
    with open(gff_path) as in_fh, open(out_path, "w") as out_fh:
        current_gene_id = None
        current_is_gene_block = False
        for line in in_fh:
            if line.startswith("#") or not line.strip():
                out_fh.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                out_fh.write(line)
                continue
            feature, attrs = fields[2], fields[8]

            if feature == "protein_coding_gene":
                m = ID_RE.search(attrs)
                gene_id = m.group(1) if m else None
                current_gene_id = gene_id
                current_is_gene_block = True
                new_id = id_mapping[gene_id]
                tag = best_tag.get(gene_id, "NA")
                attrs = add_original_id_attr(attrs, gene_id)
                attrs = replace_attr_value(attrs, "ID", gene_id, new_id)
                attrs = attrs.rstrip(";") + f";intersect={tag}"
                fields[2] = "gene"  # get_pangenes only clusters feature type "gene"
                fields[8] = attrs
                out_fh.write("\t".join(fields) + "\n")
                continue

            if feature not in CHILD_FEATURE_TYPES:
                current_is_gene_block = False
                current_gene_id = None
                out_fh.write(line)
                continue

            if current_is_gene_block:
                if feature == "pseudogenic_transcript":
                    fields[2] = "transcript"  # get_pangenes' expected child type
                if current_gene_id in id_mapping and feature in ("mRNA", "pseudogenic_transcript"):
                    attrs = replace_attr_value(attrs, "Parent", current_gene_id, id_mapping[current_gene_id])
                fields[8] = attrs
            out_fh.write("\t".join(fields) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--strains",
        required=True,
        help="Comma-separated genome names or a file with one per line - first one must be the reference, same as the rest of this pipeline's ORDER",
    )
    parser.add_argument(
        "--final-dir", required=True, help="Directory with <strain>_final_annotation.gff (this pipeline's raw output)"
    )
    parser.add_argument(
        "--out-dir", default="./tagged_annotation", help="Where to write the tagged GFFs + intersect-pair reports"
    )
    parser.add_argument(
        "--min-overlap-bp",
        type=int,
        default=20,
        help="CDS-to-CDS overlaps under this many bp are ignored entirely - boundary-rounding noise, not treated as a real overlap (default: 20)",
    )
    parser.add_argument(
        "--reference-source-tag",
        default="main",
        help="The source tag marking the reference-derived tier in this pipeline's output (default: 'main')",
    )
    parser.add_argument(
        "--id-prefix",
        default="Afu",
        help="Prefix for renamed gene IDs (default: 'Afu', this project's own convention - change for other species)",
    )
    args = parser.parse_args()

    strains = read_strains(args.strains)
    os.makedirs(args.out_dir, exist_ok=True)

    summary_rows = []
    for i, s in enumerate(strains, start=1):
        gff_path = os.path.join(args.final_dir, f"{s}_final_annotation.gff")
        out_path = os.path.join(args.out_dir, f"{s}_final_annotation_tagged.gff")
        pairs_report_path = os.path.join(args.out_dir, f"{s}_intersect_pairs.tsv")

        print(f"=== Genome {i}/{len(strains)}: {s} ===", file=sys.stderr)
        genes = parse_gff(gff_path)
        overlap_pairs = find_cds_overlaps(genes, args.min_overlap_bp)
        best_tag, partners = compute_tags(genes, overlap_pairs, args.reference_source_tag)
        id_mapping = assign_new_ids(genes, genes, i, args.id_prefix)
        write_tagged_gff(gff_path, out_path, best_tag, id_mapping)

        with open(pairs_report_path, "w") as fh:
            fh.write("gene_id\tintersect\tpartner_id\toverlap_bp\n")
            for gid, tag in best_tag.items():
                fh.writelines(f"{gid}\t{tag}\t{partner_id}\t{ov_bp}\n" for partner_id, ov_bp in partners[gid])

        n_total = len(genes)
        tag_counts = defaultdict(int)
        for gid in genes:
            tag_counts[best_tag.get(gid, "NA")] += 1
        print(
            f"  {n_total} protein-coding genes -> "
            f"NA={tag_counts['NA']}, with_main={tag_counts['with_main']}, "
            f"with_other_with_domain={tag_counts['with_other_with_domain']}, "
            f"with_other_no_domain={tag_counts['with_other_no_domain']} -> {out_path}",
            file=sys.stderr,
        )
        summary_rows.append(
            (
                s,
                n_total,
                tag_counts["NA"],
                tag_counts["with_main"],
                tag_counts["with_other_with_domain"],
                tag_counts["with_other_no_domain"],
            )
        )

    with open(os.path.join(args.out_dir, "tag_summary.tsv"), "w") as fh:
        fh.write("strain\tn_total\tn_NA\tn_with_main\tn_with_other_with_domain\tn_with_other_no_domain\n")
        fh.writelines("\t".join(str(x) for x in row) + "\n" for row in summary_rows)

    print(f"All genomes tagged (nothing removed). Outputs are in {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
