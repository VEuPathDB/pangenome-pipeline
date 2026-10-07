#!/usr/bin/env python3
"""
Post-process every strain's <genome>_final_annotation.gff (this pipeline's output) with three
filters: resolving overlapping-gene pairs, then dropping unsupported extra-copy genes, then the
rescued-coverage floor. Reads ALL genomes in --strains in one run. The overlap verdict (step 1) is
ONE global verdict per gene (by its original ID, shared across every genome that carries it) - a
gene that loses an overlap anywhere in the panel is dropped from every genome's output, even where
its own local copy never overlapped anything. The extra-copy and rescued-coverage checks (steps
2-3) judge one copy in one genome, so they drop only that copy, in that genome.

Execution order matters:

0) Panel-wide stats, from every genome's RAW genes (before any filtering): how many genomes
   carry each gene (`present_in`), how many of those are good_<genome> transfers, and - from a
   first CDS-overlap scan of every genome - how many genomes each gene has a *qualifying*
   overlap in (`collides_in`; a "minor + flank" overlap, auto-kept below, does not count).

1) Overlap resolution, decided globally per gene, using every genome's full raw gene set
   (including genes that would later fail the extra-copy/rescue-coverage checks in step 2 - an
   overlap dispute is settled before those checks ever run).

   A minor overlap (< 20% of the smaller gene's own CDS) that sits on the bigger gene's
   outermost CDS segment ("flank") is left alone: the two gene models most likely just disagree
   at the boundary, and both genes are probably real. Keep both, done.

   Everything else - minor-but-not-flank, or "substantial"/"contained" (>= 20%, sized against
   whichever gene is smaller) - needs a decision. A gene "has evidence" if it has a protein
   domain (InterPro/PFAM Dbxref hit) or a specific description (anything other than
   "hypothetical protein", "unspecified product", or empty - see GENERIC_DESCRIPTIONS below;
   --require-domain-only turns off the description half, domain-only, like this script used to
   work).

   Whichever gene has more evidence wins:
     - one gene is main, the other has no evidence -> main wins, the other is a drop candidate,
       run through the panel check below.
     - one gene is main, the other also has evidence -> keep both.
     - neither is main, only one has evidence -> that one wins, the other is a drop candidate,
       run through the panel check below.
     - neither is main, both have evidence -> keep both.
     - neither is main, neither has evidence -> no local evidence to decide by at all; no panel
       check here either - the two genes are compared directly by panel-wide good_<genome>
       transfer count instead (more wins, tied -> the longer CDS wins). This also naturally
       covers a panel-private gene: with nowhere to have transferred well into, it has 0
       good_<genome> transfers and loses to any better-supported opponent on its own.

   The panel check (cases above where one gene lost to main or to an evidence-haver) - using
   the panel-wide stats from step 0:
     - present in <= 2 genomes -> dropped (too little support to weigh against anything else).
     - present in > 2 genomes but with a qualifying overlap in only 1 of them -> NOT dropped -
       more likely a Liftoff boundary glitch specific to that one genome than real evidence the
       gene is spurious.
     - present in > 2 genomes with a qualifying overlap in >= --min-overlap-genomes (default 2)
       of them -> dropped, a real recurring conflict.

2) Extra-copy genes (extra_copy_number > 0), among whatever survived step 1, are kept only if
   they have a valid ORF - dropped otherwise, in that genome only (--extra-copy-require-valid-orf, on by default). No
   contig-size or domain check: a genuine extra copy (real tandem duplication, real gene family
   expansion) can legitimately sit on a short contig or lack a domain hit on its own.

3) Every rescued_<genome> gene, among whatever survived steps 1-2, below
   --min-rescued-coverage (default 0.5) Liftoff coverage is dropped, in that genome only. rescued_<genome>
   is this pipeline's lowest-confidence tier (round 2 of the build step - a genome's *entire*
   transfer, unfiltered, added only as a last resort) - it never had to clear a coverage/identity
   bar to get into the final annotation in the first place, unlike good_<genome>, so this is the
   one basic quality floor that tier gets.

4) Every kept protein-coding gene/pseudogene is relabeled to feature type "gene" (whatever it
   was called before), and its "pseudogenic_transcript" child (if any) to "transcript" - the one
   thing get_pangenes' own parser actually requires. Non-coding RNA genes are left exactly as
   they are.

5) Every kept gene's ID is renumbered to "<id-prefix><strain-index>g<00001, 00002, ...>". The
   original ID is kept in a new "original_ID=" attribute.

Every removed gene is logged with why, in each genome's own --out-dir/<genome>_filter_report.tsv
- nothing is silently dropped.

Known limitation - "protein domain" is only as good as the source annotation: a genome whose
ORIGINAL annotation never recorded InterPro/PFAM hits at all will have every gene sourced from
it test has_domain=False unconditionally, not because those genes lack real domains, but because
this annotation source never recorded them.

Usage:
    python3 filter_final_annotation.py \\
        --strains genome_order.txt \\
        --final-dir expanded_annotation \\
        --out-dir filtered_annotation
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.io import read_strains

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotation_common import (
    CHILD_FEATURE_TYPES,
    ID_RE,
    add_original_id_attr,
    assign_new_ids,
    classify_minor_position,
    find_cds_overlaps,
    parse_gff,
    replace_attr_value,
)

# A description this generic counts as no evidence at all - same as having no description.
GENERIC_DESCRIPTIONS = {
    "hypothetical protein",
    "hypothetical protein, conserved",
    "conserved hypothetical protein",
    "unspecified product",
    "unknown function",
    "unnamed protein product",
    "",
}


def has_evidence(gene, require_domain_only):
    if gene["has_domain"]:
        return True
    if require_domain_only:
        return False
    return gene["description"] not in GENERIC_DESCRIPTIONS


def is_qualifying_overlap(g1, g2, ov_bp):
    """A "minor + flank" overlap is auto-kept and doesn't count as a real conflict - everything
    else (minor-but-not-flank, substantial, contained) needs a decision. Returns
    (qualifies: bool, pct_of_smaller: float)."""
    smaller_bp = min(g1["cds_bp"], g2["cds_bp"])
    if smaller_bp == 0:
        return False, 0.0
    pct = ov_bp / smaller_bp
    if pct < 0.2:
        bigger, smaller = (g1, g2) if g1["cds_bp"] >= g2["cds_bp"] else (g2, g1)
        if classify_minor_position(bigger, smaller["cds_segments"]) == "flank":
            return False, pct
    return True, pct


def panel_check(gene_id, opponent_id, pct, present_in, collides_in, min_overlap_genomes, note):
    """Applied to a gene that locally lost an evidence comparison (or, in the neither-has-
    evidence case, to each side independently) before actually dropping it. Returns a reason
    string if it should be dropped, else None (survives - not enough panel-wide evidence either
    way, or an isolated collision more likely explained by Liftoff boundary noise)."""
    n_present = len(present_in.get(gene_id, ()))
    n_collide = len(collides_in.get(gene_id, ()))
    if n_present <= 2:
        return f"overlaps {opponent_id} ({round(100 * pct, 1)}% of its own CDS); {note}present in only {n_present} genome(s) panel-wide"
    if n_collide <= 1:
        return None
    if n_collide >= min_overlap_genomes:
        return f"overlaps {opponent_id} ({round(100 * pct, 1)}% of its own CDS); {note}present in {n_present} genomes, overlapping something in {n_collide} of them"
    return None


def resolve_overlaps_globally(genes_by_strain, args, present_in, n_good_transfers, collides_in, qualifying_by_strain):
    """One pass over every genome's qualifying overlap pairs, building ONE verdict per gene ID
    (by original ID, shared across the panel) - {gene_id: reason}. A gene already in verdict is
    skipped in any later pair it's part of (already decided, and verdicts never reverse)."""
    verdict = {}
    for s, genes in genes_by_strain.items():
        for gid1, gid2, ov_bp, pct in qualifying_by_strain[s]:
            if gid1 in verdict and gid2 in verdict:
                continue
            g1, g2 = genes[gid1], genes[gid2]
            g1_is_ref = g1["source"] == args.reference_source_tag
            g2_is_ref = g2["source"] == args.reference_source_tag

            if g1_is_ref and g2_is_ref:
                continue  # both reference-derived - not this pipeline's call, keep both

            if g1_is_ref or g2_is_ref:
                ref_gene, other = (g1, g2) if g1_is_ref else (g2, g1)
                if has_evidence(other, args.require_domain_only):
                    continue  # uncertain - keep both
                if other["id"] in verdict:
                    continue
                reason = panel_check(
                    other["id"],
                    ref_gene["id"],
                    pct,
                    present_in,
                    collides_in,
                    args.min_overlap_genomes,
                    f"no evidence, overlaps the reference ('{args.reference_source_tag}') gene {ref_gene['id']}; ",
                )
                if reason:
                    verdict[other["id"]] = reason
                continue

            g1_ev = has_evidence(g1, args.require_domain_only)
            g2_ev = has_evidence(g2, args.require_domain_only)

            if g1_ev != g2_ev:
                winner, loser = (g1, g2) if g1_ev else (g2, g1)
                if loser["id"] in verdict:
                    continue
                reason = panel_check(
                    loser["id"],
                    winner["id"],
                    pct,
                    present_in,
                    collides_in,
                    args.min_overlap_genomes,
                    f"no evidence, {winner['id']} has evidence; ",
                )
                if reason:
                    verdict[loser["id"]] = reason
                continue

            if g1_ev and g2_ev:
                continue  # both have evidence - uncertain, keep both

            if gid1 in verdict or gid2 in verdict:
                continue
            # neither has evidence, neither is main - no panel-presence/recurrence check here
            # (that's only for cases 3/4 above); compare the two genes directly by panel-wide
            # good_<genome> transfer count instead (tie -> longer CDS wins). This also covers
            # a panel-private gene naturally: with no genome to transfer well into, it has 0
            # good_<genome> transfers and loses to any better-supported opponent on its own.
            gt1, gt2 = n_good_transfers.get(gid1, 0), n_good_transfers.get(gid2, 0)
            if gt1 != gt2:
                winner, loser, metric = (g1, g2, f"{gt1} vs {gt2}") if gt1 > gt2 else (g2, g1, f"{gt2} vs {gt1}")
            else:
                winner, loser = (g1, g2) if g1["cds_bp"] >= g2["cds_bp"] else (g2, g1)
                metric = f"tied on good_<genome> transfers panel-wide ({gt1}), kept the longer CDS"
            verdict[loser["id"]] = (
                f"overlaps {winner['id']} ({round(100 * pct, 1)}% of its own CDS); neither has evidence, "
                f"{winner['id']} more supported panel-wide ({metric})"
            )
    return verdict


def write_filtered_gff(gff_path, out_path, genes, kept_ids, id_mapping):
    with open(gff_path) as fh:
        header_lines = [l for l in fh if l.startswith("#")]

    with open(out_path, "w") as out_fh:
        out_fh.writelines(header_lines)
        with open(gff_path) as fh:
            current_gene_kept = False
            current_gene_id = None
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                feature, attrs = fields[2], fields[8]
                if feature == "protein_coding_gene":
                    m = ID_RE.search(attrs)
                    gene_id = m.group(1) if m else None
                    current_gene_id = gene_id
                    current_gene_kept = gene_id in kept_ids
                    if current_gene_kept:
                        new_id = id_mapping[gene_id]
                        attrs = add_original_id_attr(attrs, gene_id)
                        attrs = replace_attr_value(attrs, "ID", gene_id, new_id)
                        fields[2] = "gene"  # get_pangenes only clusters feature type "gene"
                        fields[8] = attrs
                        out_fh.write("\t".join(fields) + "\n")
                    continue
                if feature not in CHILD_FEATURE_TYPES:
                    current_gene_kept = True
                    current_gene_id = None
                    out_fh.write(line)
                    continue
                if current_gene_kept:
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
        "--final-dir",
        required=True,
        help="Directory with <strain>_final_annotation.gff (this pipeline's raw output, before filtering)",
    )
    parser.add_argument("--out-dir", default="./filtered_annotation", help="Where to write the filtered GFFs + reports")
    parser.add_argument(
        "--min-overlap-bp",
        type=int,
        default=20,
        help="CDS-to-CDS overlaps under this many bp are ignored entirely - boundary-rounding noise, not treated as a real overlap (default: 20)",
    )
    parser.add_argument(
        "--min-overlap-genomes",
        type=int,
        default=2,
        help="For a gene present in >2 genomes panel-wide: it needs a qualifying overlap in at least this many of them to count as a real, recurring conflict rather than an isolated Liftoff boundary glitch (default: 2)",
    )
    parser.add_argument(
        "--require-domain-only",
        action="store_true",
        help="Only a protein domain counts as evidence - ignore description text entirely (default: off, a specific description also counts)",
    )
    parser.add_argument(
        "--reference-source-tag",
        default="main",
        help="The source tag marking the reference-derived tier in this pipeline's output (default: 'main')",
    )
    parser.add_argument(
        "--extra-copy-require-valid-orf",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Drop an extra-copy gene (extra_copy_number > 0) unless it has a valid ORF (default: on; use --no-extra-copy-require-valid-orf to keep every extra copy regardless of ORF)",
    )
    parser.add_argument(
        "--min-rescued-coverage",
        type=float,
        default=0.5,
        help="Every rescued_<genome> gene (this pipeline's lowest-confidence Liftoff tier - it never had to pass a coverage/identity bar to get in) below this Liftoff coverage is dropped outright, regardless of ORF/evidence/overlap status (default: 0.5)",
    )
    parser.add_argument(
        "--id-prefix",
        default="Afu",
        help="Prefix for renamed gene IDs (default: 'Afu', this project's own convention - change for other species)",
    )
    args = parser.parse_args()

    strains = read_strains(args.strains)
    os.makedirs(args.out_dir, exist_ok=True)

    # --- Step 0: parse every genome's raw genes, and panel-wide stats from all of them ---
    genes_by_strain = {}
    present_in = defaultdict(set)
    n_good_transfers = defaultdict(int)
    for s in strains:
        genes = parse_gff(os.path.join(args.final_dir, f"{s}_final_annotation.gff"))
        genes_by_strain[s] = genes
        for gid, g in genes.items():
            present_in[gid].add(s)
            if g["source"].startswith("good_"):
                n_good_transfers[gid] += 1
    print(
        f"Panel-wide presence built from {len(strains)} genomes: {len(present_in)} distinct gene IDs", file=sys.stderr
    )

    # First CDS-overlap scan of every genome, to know - panel-wide - how many genomes each gene
    # has a *qualifying* overlap in (minor+flank doesn't count, see module docstring).
    collides_in = defaultdict(set)
    qualifying_by_strain = {}
    for s, genes in genes_by_strain.items():
        overlap_pairs = find_cds_overlaps(genes, args.min_overlap_bp)
        qualifying = []
        for gid1, gid2, ov_bp in overlap_pairs:
            qualifies, pct = is_qualifying_overlap(genes[gid1], genes[gid2], ov_bp)
            if not qualifies:
                continue
            qualifying.append((gid1, gid2, ov_bp, pct))
            collides_in[gid1].add(s)
            collides_in[gid2].add(s)
        qualifying_by_strain[s] = qualifying
    print(
        f"Overlap scan: {sum(len(v) for v in qualifying_by_strain.values())} qualifying pairs across the panel",
        file=sys.stderr,
    )

    # --- Step 1: resolve overlaps globally - one verdict per gene ID, panel-wide ---
    verdict = resolve_overlaps_globally(
        genes_by_strain, args, present_in, n_good_transfers, collides_in, qualifying_by_strain
    )
    print(f"Overlap resolution: {len(verdict)} genes dropped panel-wide", file=sys.stderr)

    # --- Steps 2-3: extra copies, then rescued-coverage floor - only among overlap survivors ---
    # These judge one copy in one genome, so the verdict is per (genome, gene), never panel-wide.
    local_verdict = {}
    for s, genes in genes_by_strain.items():
        for gid, g in genes.items():
            if gid in verdict:
                continue
            if g["extra_copy_number"] and args.extra_copy_require_valid_orf and not g["valid_ORF"]:
                local_verdict[(s, gid)] = "extra copy with no valid ORF"
    for s, genes in genes_by_strain.items():
        for gid, g in genes.items():
            if gid in verdict or (s, gid) in local_verdict or not g["source"].startswith("rescued_"):
                continue
            if g["coverage"] is not None and g["coverage"] < args.min_rescued_coverage:
                local_verdict[(s, gid)] = (
                    f"rescued transfer, Liftoff coverage {g['coverage']:.2f} < {args.min_rescued_coverage}"
                )

    def reason_for(s, gid):
        return verdict.get(gid) or local_verdict.get((s, gid))

    # --- Steps 4-5: write every genome's output: overlap verdicts everywhere, steps 2-3 per genome ---
    summary_rows = []
    for i, s in enumerate(strains, start=1):
        genes = genes_by_strain[s]
        gff_path = os.path.join(args.final_dir, f"{s}_final_annotation.gff")
        out_path = os.path.join(args.out_dir, f"{s}_final_annotation_filtered.gff")
        report_path = os.path.join(args.out_dir, f"{s}_filter_report.tsv")

        kept_ids = {gid for gid in genes if not reason_for(s, gid)}
        id_mapping = assign_new_ids(genes, kept_ids, i, args.id_prefix)
        write_filtered_gff(gff_path, out_path, genes, kept_ids, id_mapping)

        with open(report_path, "w") as rep_fh:
            rep_fh.write("gene_id\tsource\tseqname\treason\n")
            for gid, g in genes.items():
                if reason_for(s, gid):
                    rep_fh.write(f"{gid}\t{g['source']}\t{g['seqname']}\t{reason_for(s, gid)}\n")

        n_total = len(genes)
        n_removed = sum(1 for gid in genes if reason_for(s, gid))
        print(
            f"  {s}: {n_total} protein-coding genes -> {n_total - n_removed} kept, "
            f"{n_removed} removed ({round(100 * n_removed / n_total, 1) if n_total else 0}%) -> {out_path}",
            file=sys.stderr,
        )
        summary_rows.append((s, n_total, n_total - n_removed, n_removed))

    with open(os.path.join(args.out_dir, "filter_summary.tsv"), "w") as fh:
        fh.write("strain\tn_original\tn_kept\tn_removed\n")
        for s, n_total, n_kept, n_removed in summary_rows:
            fh.write(f"{s}\t{n_total}\t{n_kept}\t{n_removed}\n")

    print(f"All genomes filtered. Outputs are in {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
