#!/usr/bin/env python3
"""
Intersect curated TE annotations (RepeatMasker) with gene models (GFF3) to
see how much of each gene's transcribed sequence overlaps a TE - one signal
for flagging gene calls that are likely mis-annotated TEs (or genuinely
TE-domesticated genes).

For each genome this script:
  1. Parses the RepeatMasker .out file the same way as summarize_te_content.py
     (dropping Simple_repeat/Low_complexity hits, merging overlapping hits
     per contig into non-redundant TE intervals - both overall and per TE
     class, so a gene's dominant overlapping TE class can be reported).
  2. Parses the matching GFF3. For protein_coding genes, the "transcript"
     span is the union of CDS intervals across all of the gene's mRNAs.
     Pseudogenes are handled separately: CDS union if the pseudogene has one
     (some do - a degenerate/frameshifted copy), else exon union, else the
     gene's own genomic span - because a pseudogene not uncommonly has no
     CDS feature at all, unlike a protein_coding gene.

     tRNA / rRNA genes are not specially computed here: TE overlap on a
     ~70-90 bp tRNA is rarely a meaningful signal the way it is for a
     protein-coding CDS or a pseudogene, so they are just counted in the
     skipped-biotype warning like any other non-target biotype.
  3. Intersects each gene's transcript-interval union with the TE intervals
     on the same contig -> TE_overlap_bp / TE_overlap_pct, plus which TE
     class contributes the most overlap.

Outputs one combined per-gene TSV across all genomes, threshold-sensitivity
summaries, and exploratory plots - kept separate for protein_coding genes
(the main CDS-vs-TE question, and where the gene-model-check panels apply)
vs. pseudogenes (a different span definition, no gene-model-check panels).

Usage:
    python3 te_gene_overlap.py --te-indir repeatmasker_curated_annotation_010926 \
        --gff-indir gff_originals --outdir te_annotation_integration

Expected input layout:
  --te-indir : one R_<genome_id>/ directory per genome, each containing
               exactly one RepeatMasker .out file (same convention as
               summarize_te_content.py --indir; genome_id can be anything,
               e.g. an accession alone or accession+assembly-name).
  --gff-indir: one GFF3 file per genome (any *.gff / *.gff3 filename).

No particular accession format (e.g. NCBI "GCA_...") is assumed anywhere.
Each TE genome_id is matched to a GFF file by, in order: (1) an entry in
--genome-map if supplied, (2) an exact filename-stem match, (3) a unique
prefix match in either direction (handles the common case where one side's
name is the other's name plus a suffix, e.g. genome_id "X_assemblyname" vs.
GFF file "X.gff"). A genome_id that still can't be matched unambiguously is
skipped with a warning - use --genome-map (a 2-column TSV of
genome_id<TAB>gff_path) to resolve those explicitly.
"""

import argparse
import bisect
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import sys

import matplotlib.pyplot as plt
import pandas as pd
from summarize_te_content import (
    CATEGORICAL,
    INK_PRIMARY,
    INK_SECONDARY,
    REMOVE_CLASSES,
    SURFACE,
    find_genome_dirs,
    parse_out,
    style_axes,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.gff import (  # GENE_FTYPE_TO_BIOTYPE, parse_attributes: re-exported to sibling scripts
    GENE_FTYPE_TO_BIOTYPE,
    gene_description,
    gff_seqids,
    infer_biotype,
    parse_attributes,
)

THRESHOLDS = [1, 5, 10, 20, 50, 80, 90, 100]


def index_gff_dir(gff_indir: Path) -> dict:
    """filename stem (no extension) -> path, for every .gff/.gff3 file."""
    return {p.stem: p for p in sorted(gff_indir.glob("*.gff")) + sorted(gff_indir.glob("*.gff3"))}


def load_genome_map(path) -> dict:
    """Optional 2-column TSV (genome_id<TAB>gff_path) for explicit
    TE-genome-to-GFF matches, for naming conventions the automatic matching
    below can't resolve on its own."""
    mapping = {}
    with Path(path).open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genome_id, gff_path = line.split("\t")
            mapping[genome_id] = Path(gff_path)
    return mapping


ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")


def find_gff_for_genome(
    genome_id: str, gff_by_stem: dict, genome_map: dict, te_contigs: set | None = None, seqid_cache: dict | None = None
):
    """Match a RepeatMasker genome_id (its R_<genome_id>/ directory name) to a
    GFF file, without assuming any naming convention. Tried in order:
    (1) an explicit --genome-map entry;
    (2) an exact filename-stem match, case-insensitive;
    (3) a unique prefix match in either direction, case-insensitive (genome_id
        "X_assemblyname" vs GFF "X.gff");
    (4) a GC[AF]_ accession that appears in both names;
    (5) which GFF's seqids (column 1) overlap the RepeatMasker .out contigs -
        naming-independent, this catches the case where the two sides share no
        name at all (e.g. TE folder "RM_AF293", GFF "FungiDB-68_Af293.gff").
    Returns None if nothing resolves."""
    if genome_id in genome_map:
        return genome_map[genome_id]
    if genome_id in gff_by_stem:
        return gff_by_stem[genome_id]

    lower = {stem.lower(): stem for stem in gff_by_stem}
    gl = genome_id.lower()
    if gl in lower:
        return gff_by_stem[lower[gl]]

    candidates = [stem for stem in gff_by_stem if gl.startswith(stem.lower()) or stem.lower().startswith(gl)]
    if len(candidates) == 1:
        return gff_by_stem[candidates[0]]
    if len(candidates) > 1:
        candidates.sort(key=len, reverse=True)
        if len(candidates[0]) > len(candidates[1]):
            return gff_by_stem[candidates[0]]

    m = ACC_RE.search(genome_id)
    if m:
        acc = m.group(0)
        acc_hits = [stem for stem in gff_by_stem if acc in stem]
        if len(acc_hits) == 1:
            return gff_by_stem[acc_hits[0]]

    if te_contigs:
        cache = seqid_cache if seqid_cache is not None else {}
        best_path, best_frac = None, 0.0
        for stem, path in gff_by_stem.items():
            sids = cache.get(stem)
            if sids is None:
                sids = cache[stem] = gff_seqids(path)
            if not sids:
                continue
            frac = len(te_contigs & sids) / len(te_contigs)
            if frac > best_frac:
                best_path, best_frac = path, frac
        if best_frac >= 0.5:
            return best_path

    return None


def merge_intervals(intervals):
    """Sort + merge overlapping/touching (start,end) 1-based inclusive
    intervals. For a merged, sorted-by-start, non-overlapping list, both
    starts and ends end up strictly increasing (no interval can be nested
    inside another - they'd have been merged together) - MergedTrack below
    relies on that to binary-search instead of scanning."""
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = []
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


class MergedTrack:
    """A merged, non-overlapping interval set on one contig, indexed for
    fast overlap-bp queries against a (typically much smaller) query
    interval list."""

    def __init__(self, merged_intervals):
        self.intervals = merged_intervals
        self.starts = [s for s, _ in merged_intervals]
        self.ends = [e for _, e in merged_intervals]

    def overlap_bp(self, query_intervals) -> int:
        if not self.intervals or not query_intervals:
            return 0
        total = 0
        for s, e in query_intervals:
            lo = bisect.bisect_left(self.ends, s)
            hi = bisect.bisect_right(self.starts, e) - 1
            for idx in range(lo, hi + 1):
                bs, be = self.intervals[idx]
                ov_s, ov_e = max(s, bs), min(e, be)
                if ov_e >= ov_s:
                    total += ov_e - ov_s + 1
        return total

    def overlap_ranges(self, s: int, e: int) -> list:
        """Same matching as overlap_bp, but for a single query interval and
        returning the actual overlapping (start,end) sub-ranges instead of
        just their summed length - used to locate *where* within a gene an
        overlap falls, not just how much of it overlaps."""
        if not self.intervals:
            return []
        lo = bisect.bisect_left(self.ends, s)
        hi = bisect.bisect_right(self.starts, e) - 1
        ranges = []
        for idx in range(lo, hi + 1):
            bs, be = self.intervals[idx]
            ov_s, ov_e = max(s, bs), min(e, be)
            if ov_e >= ov_s:
                ranges.append((ov_s, ov_e))
        return ranges


def build_te_tracks(out_path: Path):
    """Per contig: a MergedTrack of all kept TE hits, a dict of MergedTrack
    per TE class/family (RepeatMasker's class/family column), a dict of
    MergedTrack per TE name (RepeatMasker's matching-repeat column, e.g.
    "GYMAG1_Pory") - used to report which TE class and which specific named
    TE model dominate a gene's overlap - and a list of de-fragmented element
    spans (see build_element_spans) - used to report how large the actual
    TE copy behind an overlap is, not just how many bp of it overlap."""
    hits = parse_out(out_path)
    kept = [h for h in hits if h[4].lower() not in REMOVE_CLASSES]

    by_contig = defaultdict(list)
    by_contig_class: defaultdict[str, defaultdict[str, list]] = defaultdict(lambda: defaultdict(list))
    by_contig_name: defaultdict[str, defaultdict[str, list]] = defaultdict(lambda: defaultdict(list))
    for seq, start, end, name, cls, seq_left, pct_div, rm_id in kept:
        by_contig[seq].append((start, end))
        by_contig_class[seq][cls].append((start, end))
        by_contig_name[seq][name].append((start, end))

    te_all = {seq: MergedTrack(merge_intervals(ivs)) for seq, ivs in by_contig.items()}
    te_by_class = {
        seq: {cls: MergedTrack(merge_intervals(ivs)) for cls, ivs in classes.items()}
        for seq, classes in by_contig_class.items()
    }
    te_by_name = {
        seq: {name: MergedTrack(merge_intervals(ivs)) for name, ivs in names.items()}
        for seq, names in by_contig_name.items()
    }
    element_spans = build_element_spans(kept)
    return te_all, te_by_class, te_by_name, element_spans


def build_element_spans(kept_hits):
    """Per contig: a start-sorted list of (start, end, cls, name, pct_div)
    for each de-fragmented TE *element* (a full inserted copy), not each
    RepeatMasker hit line. RepeatMasker splits one element into several hit
    lines when a later, unrelated insertion lands inside it (its own ID
    column groups those pieces back together on a contig); this takes each
    such group's full span (first piece's start to last piece's end) and its
    length-weighted average % divergence from consensus (a multi-fragment
    element's pieces can each score slightly differently).

    This is the same defragmentation summarize_te_content.py's
    num_TE_fragments count relies on, resurfaced here for a different use:
    reporting how large the actual TE copy behind a gene's overlap is - a
    gene barely touched (low overlap bp) by a 2 kb element is a much
    stronger "this is a real TE" signal than being barely touched by a
    ~50 bp one, which could just as easily be a spurious/degenerate
    fragment. Unlike merge_intervals's output, these spans are NOT
    flattened into non-overlapping blocks: a nested insertion's span sits
    entirely inside its host element's, by construction, and both need to
    stay visible as separate elements."""
    raw = defaultdict(dict)  # seq -> rm_id -> [lo, hi, cls, name, div_num, div_den]
    for seq, start, end, name, cls, seq_left, pct_div, rm_id in kept_hits:
        frag_len = end - start + 1
        rec = raw[seq].get(rm_id)
        if rec is None:
            raw[seq][rm_id] = [start, end, cls, name, pct_div * frag_len, frag_len]
        else:
            rec[0] = min(rec[0], start)
            rec[1] = max(rec[1], end)
            rec[4] += pct_div * frag_len
            rec[5] += frag_len
    return {
        seq: sorted((lo, hi, cls, name, div_num / div_den) for lo, hi, cls, name, div_num, div_den in spans.values())
        for seq, spans in raw.items()
    }


def overlapping_elements(element_list, query_intervals):
    """element_list: this contig's build_element_spans() output - start-
    sorted, entries may be nested/overlapping each other. query_intervals:
    a gene's merged CDS blocks. Returns every element that overlaps at
    least one query block, as
    (start, end, cls, name, overlap_bp_with_gene, pct_div) - length =
    end - start + 1 is that element's full (de-fragmented) size, regardless
    of how much of it falls inside the gene."""
    results = []
    for s, e in query_intervals:
        for el_s, el_e, cls, name, div in element_list:
            if el_s > e:
                break  # start-sorted: no later element can start <= e either
            if el_e < s:
                continue
            ov = min(e, el_e) - max(s, el_s) + 1
            if ov > 0:
                results.append((el_s, el_e, cls, name, ov, div))
    return results


def parse_gff(gff_path: Path):
    """Single pass: gene/pseudogene -> metadata (incl. gene_biotype and, when
    present, a free-text description), RNA feature -> parent gene (plus its
    own description, for the common case where only the mRNA/CDS carries
    one), and CDS/exon intervals keyed by their parent RNA (resolved to
    genes afterwards, so this doesn't depend on parents appearing before
    children in the file)."""
    gene_meta = {}
    rna_to_gene = {}
    rna_description = {}
    cds_by_rna = defaultdict(list)
    exon_by_rna = defaultdict(list)

    with gff_path.open() as fh:
        for line in fh:
            if line[0] == "#" or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _source, ftype, start, end, _score, strand, _phase, attr_str = fields

            if ftype in GENE_FTYPE_TO_BIOTYPE:
                attrs = parse_attributes(attr_str)
                gid = attrs.get("ID")
                if gid is None:
                    continue
                gene_meta[gid] = {
                    "contig": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_biotype": infer_biotype(ftype, attrs),
                    "description": gene_description(attrs),
                }
            elif ftype in ("mRNA", "tRNA", "rRNA"):
                attrs = parse_attributes(attr_str)
                rid = attrs.get("ID")
                parent = attrs.get("Parent")
                if rid and parent:
                    rna_to_gene[rid] = parent
                if rid:
                    desc = gene_description(attrs)
                    if desc:
                        rna_description[rid] = desc
            elif ftype == "CDS":
                parent = parse_attributes(attr_str).get("Parent")
                if parent:
                    cds_by_rna[parent].append((int(start), int(end)))
            elif ftype == "exon":
                parent = parse_attributes(attr_str).get("Parent")
                if parent:
                    exon_by_rna[parent].append((int(start), int(end)))

    return gene_meta, rna_to_gene, rna_description, cds_by_rna, exon_by_rna


def aggregate_by_gene(rna_to_gene: dict, feature_by_rna: dict) -> dict:
    by_gene = defaultdict(list)
    for rna_id, intervals in feature_by_rna.items():
        gene_id = rna_to_gene.get(rna_id)
        if gene_id:
            by_gene[gene_id].extend(intervals)
    return by_gene


def strip_prefix(feature_id: str, prefix: str) -> str:
    return feature_id.removeprefix(prefix)


MIN_OVERLAP_BP = 50  # a raw TE overlap below this many bp is treated as no overlap at all
EDGE_FRAC = 0.15  # first/last 15% of the spliced transcript counts as an "edge" zone
MIN_LEN_FOR_POSITION = 150  # bp; below this, 5'/3'/internal isn't a meaningful distinction
EDGE_DOMINANCE_FRAC = 0.7  # >=70% of a gene's overlap bp sitting in one zone -> that zone "dominates"


def classify_overlap_position(merged: list, strand: str, track, overlap_bp: int, transcript_bp: int):
    """Where within a gene's spliced (5'->3') transcript its TE overlap
    falls - not just how much overlaps, but whether that overlap sits at
    one end (a candidate sign of a chimeric gene model: the caller fused a
    neighboring TE onto the true gene boundary) versus spread through the
    middle (more consistent with the annotated CDS itself being TE-derived,
    or a coincidental partial match).

    merged: the gene's own merged CDS blocks, genomic coords, sorted
    ascending (as merge_intervals returns them). track: that contig's
    MergedTrack of TE calls (te_all), already known to have overlap_bp > 0
    with this gene.

    Returns (position_label, overlap_5prime_bp, overlap_3prime_bp,
    overlap_internal_bp). position_label is one of:
      "spans_gene"      - already >=80% overlap; not a partial/edge case
      "short_transcript" - too short (< MIN_LEN_FOR_POSITION) to classify
      "5_prime" / "3_prime" - >=70% of the overlap sits in that terminal zone
      "internal"         - >=70% of the overlap sits away from both edges
      "mixed"            - no single zone dominates
    """
    if overlap_bp / transcript_bp >= 0.80:
        return "spans_gene", None, None, None
    if transcript_bp < MIN_LEN_FOR_POSITION:
        return "short_transcript", None, None, None

    # 5'->3' block order: ascending genomic order for +, reversed for - (a
    # minus-strand gene's 5' end is at the highest genomic coordinate) -
    # and, symmetrically, position-within-a-block also runs backwards.
    blocks = merged if strand != "-" else list(reversed(merged))
    ov5_bp = ov3_bp = 0
    edge_bp = max(1, round(EDGE_FRAC * transcript_bp))
    three_prime_zone_start = transcript_bp - edge_bp
    offset = 0
    for g_s, g_e in blocks:
        length = g_e - g_s + 1
        for ov_s, ov_e in track.overlap_ranges(g_s, g_e):
            if strand == "-":
                sp_s, sp_e = offset + (g_e - ov_e), offset + (g_e - ov_s)
            else:
                sp_s, sp_e = offset + (ov_s - g_s), offset + (ov_e - g_s)
            ov5_bp += max(0, min(sp_e, edge_bp - 1) - sp_s + 1) if sp_s < edge_bp else 0
            ov3_bp += max(0, sp_e - max(sp_s, three_prime_zone_start) + 1) if sp_e >= three_prime_zone_start else 0
        offset += length

    ov_internal_bp = overlap_bp - ov5_bp - ov3_bp
    if ov5_bp >= EDGE_DOMINANCE_FRAC * overlap_bp:
        label = "5_prime"
    elif ov3_bp >= EDGE_DOMINANCE_FRAC * overlap_bp:
        label = "3_prime"
    elif ov_internal_bp >= EDGE_DOMINANCE_FRAC * overlap_bp:
        label = "internal"
    else:
        label = "mixed"
    return label, ov5_bp, ov3_bp, ov_internal_bp


def compute_gene_rows(
    genome_id: str,
    genome_label: str,
    gene_meta,
    rna_to_gene,
    rna_description,
    cds_by_rna,
    exon_by_rna,
    te_all,
    te_by_class,
    te_by_name,
    element_spans,
    min_overlap_bp=MIN_OVERLAP_BP,
):
    """Returns (gene_rows, transcript_rows, skipped_biotypes). gene_rows is
    one row per gene: protein_coding genes use the CDS union as their span;
    pseudogenes use CDS union if present, else exon union, else the gene's
    own genomic span (a pseudogene not uncommonly has no CDS feature at
    all). transcript_rows drills into the per-isoform CDS of protein_coding
    genes that (a) have more than one mRNA and (b) have TE_overlap_bp > 0 at
    the gene level - i.e. genes where it's worth checking whether every
    isoform is TE-overlapping, or only some (which would point at a
    specific transcript model being the spurious/TE-derived one rather than
    the gene as a whole).
    skipped_biotypes is a Counter of gene_biotype -> count for genes whose
    gene_biotype was neither "protein_coding" nor "pseudogene" (this
    includes tRNA/rRNA genes - TE overlap on those is not computed) - these
    genes are silently absent from gene_rows, so the caller can warn instead
    of that going unnoticed (this matters most on a GFF from a different
    source/pipeline than this project's, which may not spell
    "protein_coding" the same way NCBI's annotwriter does).

    genome_id is the full RepeatMasker-side id (kept in the TSVs, consistent
    with the rest of the pipeline's per-genome tables); genome_label is the
    matched GFF file's stem, used only for compact plot display - it's
    whatever the GFF happened to be named, not a parsed/assumed accession.

    min_overlap_bp: any gene whose raw TE_overlap_bp comes out below this is
    treated as having 0 overlap - a real base coincidentally touching a TE
    (e.g. right at a contig/gene boundary) isn't evidence of anything on
    its own, no matter how large that TE element happens to be elsewhere."""
    cds_by_gene = aggregate_by_gene(rna_to_gene, cds_by_rna)
    exon_by_gene = aggregate_by_gene(rna_to_gene, exon_by_rna)
    gene_to_rnas = defaultdict(list)
    for rna_id, gene_id in rna_to_gene.items():
        gene_to_rnas[gene_id].append(rna_id)

    rows = []
    transcript_rows = []
    skipped_biotypes: Counter[str] = Counter()
    for gene_id, meta in gene_meta.items():
        biotype = meta["gene_biotype"]
        if biotype == "protein_coding":
            category = "protein_coding"
            raw_intervals = cds_by_gene.get(gene_id)
        elif biotype == "pseudogene":
            category = "pseudogene"
            # CDS if this pseudogene has one (a degenerate/frameshifted copy
            # sometimes still gets one), else its exons, else - if it has no
            # transcript-level features at all, common for a bare pseudogene
            # record - fall back to the gene's own genomic span.
            raw_intervals = cds_by_gene.get(gene_id) or exon_by_gene.get(gene_id) or [(meta["start"], meta["end"])]
        else:
            skipped_biotypes[biotype] += 1
            continue
        if not raw_intervals:
            continue

        merged = merge_intervals(raw_intervals)
        transcript_bp = sum(e - s + 1 for s, e in merged)
        if transcript_bp == 0:
            continue

        contig = meta["contig"]
        track = te_all.get(contig)
        overlap_bp = track.overlap_bp(merged) if track else 0
        if overlap_bp < min_overlap_bp:
            overlap_bp = 0
        overlap_pct = 100.0 * overlap_bp / transcript_bp

        dominant_class, dominant_class_bp = "", 0
        dominant_name, dominant_name_bp = "", 0
        overlap_position = ""
        max_fragment_bp, max_fragment_name, max_fragment_overlap_bp = 0, "", 0
        dominant_overlap_bp, dominant_overlap_div = 0, None
        n_te_fragments = 0
        if overlap_bp > 0:
            for cls, ctrack in te_by_class.get(contig, {}).items():
                bp = ctrack.overlap_bp(merged)
                if bp > dominant_class_bp:
                    dominant_class_bp, dominant_class = bp, cls
            for name, ntrack in te_by_name.get(contig, {}).items():
                bp = ntrack.overlap_bp(merged)
                if bp > dominant_name_bp:
                    dominant_name_bp, dominant_name = bp, name
            if category == "protein_coding":
                overlap_position, _ov5_bp, _ov3_bp, _ov_internal_bp = classify_overlap_position(
                    merged, meta["strand"], track, overlap_bp, transcript_bp
                )

            elements = overlapping_elements(element_spans.get(contig, []), merged)
            n_te_fragments = len(elements)
            for el_s, el_e, el_cls, el_name, el_overlap_bp, el_div in elements:
                el_len = el_e - el_s + 1
                if el_len > max_fragment_bp:
                    max_fragment_bp = el_len
                    max_fragment_name = el_name
                    max_fragment_overlap_bp = el_overlap_bp
                # dominant_overlap_div: % divergence of the element (a
                # de-fragmented, rm_id-grouped TE copy - see
                # build_element_spans) that actually explains the most of
                # this gene's overlap, not just the largest element that
                # happens to touch it - so a gene barely clipped by a huge,
                # near-identical element doesn't get that element's low
                # divergence misreported as "this gene = pristine TE copy".
                if el_overlap_bp > dominant_overlap_bp:
                    dominant_overlap_bp = el_overlap_bp
                    dominant_overlap_div = el_div

        gene_rnas = gene_to_rnas.get(gene_id, [])
        n_transcripts = len(gene_rnas)
        # the gene feature's own description if it has one (VEuPathDB-style);
        # else the first of its transcripts that has one (NCBI puts "product"
        # on the mRNA, not the gene)
        description = meta["description"]
        if not description:
            for rna_id in gene_rnas:
                if rna_id in rna_description:
                    description = rna_description[rna_id]
                    break
        rows.append(
            {
                "genome_id": genome_id,
                "genome_label": genome_label,
                "gene_id": strip_prefix(gene_id, "gene-"),
                "gene_description": description,
                "gene_category": category,
                "gene_biotype": biotype,
                "contig": contig,
                "start": meta["start"],
                "end": meta["end"],
                "strand": meta["strand"],
                "n_transcripts": n_transcripts,
                "transcript_bp": transcript_bp,
                "TE_overlap_bp": overlap_bp,
                "TE_overlap_pct": overlap_pct,
                "TE_divergence_pct": dominant_overlap_div,
                "dominant_TE_class": dominant_class,
                "dominant_TE_class_pct_of_overlap": (100.0 * dominant_class_bp / overlap_bp if overlap_bp else 0.0),
                "dominant_TE_name": dominant_name,
                "dominant_TE_name_pct_of_overlap": (100.0 * dominant_name_bp / overlap_bp if overlap_bp else 0.0),
                "TE_overlap_position": overlap_position,
                "check_gene_model": overlap_position in ("5_prime", "3_prime"),
                "n_TE_fragments_overlapping": n_te_fragments,
                "TE_fragment_max_bp": max_fragment_bp,
                "TE_fragment_max_bp_name": max_fragment_name,
                "TE_fragment_max_bp_overlap_bp": max_fragment_overlap_bp,
            }
        )

        if category == "protein_coding" and n_transcripts > 1 and overlap_bp > 0:
            for rna_id in gene_to_rnas[gene_id]:
                tx_intervals = cds_by_rna.get(rna_id)
                if not tx_intervals:
                    continue
                tx_merged = merge_intervals(tx_intervals)
                tx_bp = sum(e - s + 1 for s, e in tx_merged)
                if tx_bp == 0:
                    continue
                tx_overlap_bp = track.overlap_bp(tx_merged) if track else 0
                transcript_rows.append(
                    {
                        "genome_id": genome_id,
                        "genome_label": genome_label,
                        "gene_id": strip_prefix(gene_id, "gene-"),
                        "transcript_id": strip_prefix(rna_id, "rna-"),
                        "n_transcripts": n_transcripts,
                        "gene_TE_overlap_pct": overlap_pct,
                        "transcript_bp": tx_bp,
                        "transcript_TE_overlap_bp": tx_overlap_bp,
                        "transcript_TE_overlap_pct": 100.0 * tx_overlap_bp / tx_bp,
                    }
                )

    return rows, transcript_rows, skipped_biotypes


def threshold_table(df: pd.DataFrame, thresholds=THRESHOLDS) -> pd.DataFrame:
    """Cumulative >=X% counts (as before) plus one non-cumulative "20-80%"
    band column - the exact count of genes in the proposed "inspect" zone
    (20% <= overlap < 80%), which the cumulative columns alone don't give
    directly (>=20% also includes every gene that's >=80%)."""

    def band_20_80(sub):
        return int(((sub["TE_overlap_pct"] >= 20) & (sub["TE_overlap_pct"] < 80)).sum())

    rows = []
    for genome_id, sub in df.groupby("genome_id"):
        row = {"genome_id": genome_id, "genome_label": sub["genome_label"].iloc[0], "n_genes": len(sub)}
        for t in thresholds:
            row[f">={t}%"] = int((sub["TE_overlap_pct"] >= t).sum())
        row["20-80%"] = band_20_80(sub)
        rows.append(row)
    total_row = {"genome_id": "ALL", "genome_label": "ALL", "n_genes": len(df)}
    for t in thresholds:
        total_row[f">={t}%"] = int((df["TE_overlap_pct"] >= t).sum())
    total_row["20-80%"] = band_20_80(df)
    rows.append(total_row)
    return pd.DataFrame(rows)


def annotate_n(ax, n):
    """Small 'n = <count>' label centered at the top of the plot, so a
    reader never has to guess how many genes/points a plot summarizes."""
    ax.text(0.5, 0.95, f"n = {n:,}", transform=ax.transAxes, ha="center", va="top", fontsize=9, color=INK_SECONDARY)


def plot_hist(values, title, dest: Path, bins=50, log_y=False, color=CATEGORICAL[0]):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(values, bins=bins, range=(0, 100), color=color, zorder=3)
    style_axes(ax)
    if log_y:
        ax.set_yscale("log")
        ax.set_ylabel("Number of genes (log scale)")
    else:
        ax.set_ylabel("Number of genes")
    ax.set_xlabel("TE overlap (% of transcript)")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    annotate_n(ax, len(values))
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_hist_by_genome(df: pd.DataFrame, title: str, dest: Path, bins=20, log_y=True, ncols=3):
    genomes = sorted(df["genome_label"].unique())
    nrows = -(-len(genomes) // ncols)
    # sharey=True forces every panel onto the same y-axis limits, so bar
    # heights are directly comparable across genomes instead of each panel
    # auto-scaling to its own max (which would visually flatten a genome
    # with few TE-overlapping genes to look the same as one with many).
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 3.2 * nrows), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, genome_label in zip(axes, genomes):
        vals = df.loc[df["genome_label"] == genome_label, "TE_overlap_pct"]
        counts, bin_edges, patches = ax.hist(vals, bins=bins, range=(0, 100), color=CATEGORICAL[0], zorder=3)
        if log_y:
            ax.set_yscale("log")
        style_axes(ax)
        ax.set_title(f"{genome_label}\nn={len(vals):,}", fontsize=11, color=INK_SECONDARY, pad=6)
        ax.tick_params(labelsize=9)
        # only label bars in the high-overlap half (>=50%) - that's the
        # "probably a TE" region and the one worth reading exact counts
        # off of; labeling every small low-overlap bar too just added
        # clutter without answering a different question.
        for count, edge, patch in zip(counts, bin_edges[:-1], patches):
            if count > 0 and edge >= 50:
                # offset in points (not data units) so the gap above the bar
                # stays constant regardless of the panel's y-scale
                ax.annotate(
                    f"{int(count)}",
                    xy=(patch.get_x() + patch.get_width() / 2, count),
                    xytext=(0, 5),
                    textcoords="offset points",
                    rotation=45,
                    ha="left",
                    va="bottom",
                    fontsize=8,
                    color=INK_SECONDARY,
                )
    # headroom for the angled count labels above the tallest bar, applied
    # once since sharey means every panel uses the same limits anyway
    top = axes[0].get_ylim()[1]
    axes[0].set_ylim(top=top * (2.2 if log_y else 1.25))
    for ax in axes[len(genomes) :]:
        ax.axis("off")
    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_threshold_heatmap(thresh_df: pd.DataFrame, title: str, dest: Path):
    """Per-genome breakdown of the pooled threshold bar chart above: rows =
    genome, columns = TE-overlap threshold, cell = gene count at/above that
    threshold in that genome. Sorted by the >=80% column so the genomes
    with the heaviest TE-gene burden end up at the top."""
    df = thresh_df[thresh_df["genome_label"] != "ALL"].copy()
    cols = [c for c in df.columns if c.startswith(">=")]
    df = df.sort_values(">=80%", ascending=True)
    labels = df["genome_label"].tolist()
    matrix = df[cols].to_numpy(dtype=float)
    # ">=100%" is a misleading label - 100% is the maximum possible overlap,
    # so this column is really just "=100%"; every other column keeps its
    # ">=" cumulative label since those genuinely are lower bounds.
    col_labels = ["100%" if c == ">=100%" else c for c in cols]

    fig, ax = plt.subplots(figsize=(8, max(4, 0.32 * len(df) + 1.5)))
    im = ax.imshow(matrix, cmap="Blues", aspect="auto", origin="lower")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    vmax = matrix.max() if matrix.size else 1
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            color = SURFACE if val > vmax * 0.6 else INK_PRIMARY
            ax.text(j, i, f"{int(val)}", ha="center", va="center", fontsize=7, color=color)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("Number of genes", color=INK_SECONDARY)
    cb.ax.tick_params(colors=INK_SECONDARY, labelsize=8)
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def plot_band_bar(thresh_df: pd.DataFrame, column: str, xlabel: str, title: str, dest: Path, color=CATEGORICAL[3]):
    """One horizontal bar per genome, for any single column of
    threshold_table()'s output - e.g. the non-cumulative "20-80%" band, or
    one of the cumulative ">=X%" columns (like ">=80%", the eliminated
    band) - showing how many of that genome's genes fall in it."""
    df = thresh_df[thresh_df["genome_label"] != "ALL"].copy()
    values = df.set_index("genome_label")[column].sort_values()
    fig, ax = plt.subplots(figsize=(7, 0.32 * len(values) + 1.5))
    ax.barh(values.index, values.to_numpy(), color=color, height=0.65, zorder=3)
    style_axes(ax, horizontal=True)
    for y, v in enumerate(values.values):
        ax.annotate(
            f"{int(v)}",
            xy=(v, y),
            xytext=(5, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=8,
            color=INK_SECONDARY,
        )
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.text(
        0.98,
        0.03,
        f"total = {int(values.sum()):,}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color=INK_SECONDARY,
    )
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


OVERLAP_POSITION_LABELS = {
    "5_prime": "5_prime (edge, near start)",
    "3_prime": "3_prime (edge, near end)",
    "internal": "internal (away from both ends)",
    "mixed": "mixed (no single zone ≥ 70% of the overlap)",
    "short_transcript": f"short_transcript (<{MIN_LEN_FOR_POSITION} bp, not classified)",
}
OVERLAP_POSITION_COLORS = {
    "5_prime": CATEGORICAL[1],
    "3_prime": CATEGORICAL[7],
    "internal": CATEGORICAL[0],
    "mixed": CATEGORICAL[3],
    "short_transcript": "#c3c2b7",
}


def plot_overlap_position_bar(pc_df: pd.DataFrame, title: str, dest: Path):
    """Bar chart of TE_overlap_position categories for protein-coding genes
    in the 20-80% "inspect" band specifically (same scope as the rest of
    this section's threshold discussion - genes at >=80% - "spans_gene" -
    are already handled by the main threshold, and genes below 20% aren't
    the "inspect" set either). 5_prime/3_prime are highlighted since those
    are the check_gene_model candidates - overlap concentrated at one end
    of the gene, a candidate sign the gene model itself is chimeric (fused
    onto a neighboring TE) rather than the CDS being TE-derived."""
    df = pc_df[(pc_df["TE_overlap_pct"] >= 20) & (pc_df["TE_overlap_pct"] < 80)]
    order = list(OVERLAP_POSITION_COLORS)
    counts = df["TE_overlap_position"].value_counts().reindex(order, fill_value=0)
    fig, ax = plt.subplots(figsize=(10.5, 5))
    for x, cat in enumerate(order):
        v = counts[cat]
        ax.bar(x, v, color=OVERLAP_POSITION_COLORS[cat], zorder=3, label=OVERLAP_POSITION_LABELS[cat])
        ax.annotate(
            f"{int(v)}",
            xy=(x, v),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            color=INK_SECONDARY,
        )
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    style_axes(ax)
    ax.set_ylabel("Number of genes")
    ax.set_xlabel("Where the TE overlap falls along the gene (20-80% 'inspect' band only)")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8, labelcolor=INK_SECONDARY)
    # top-center (annotate_n's usual spot) lands on the tallest bars here -
    # top-left instead, above the short 5_prime bar, stays clear
    ax.text(
        0.02,
        0.95,
        f"n = {int(counts.sum()):,}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color=INK_SECONDARY,
    )
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


OVERLAP_PCT_BINS = [(0, 5), (5, 10), (10, 15), (15, 20), (20, 30), (30, 40), (40, 50), (50, 65), (65, 80)]


def plot_overlap_position_by_pct(pc_df: pd.DataFrame, title: str, dest: Path):
    """Same TE_overlap_position categories as plot_overlap_position_bar (any
    0-80% overlap, not just the 20-80% "inspect" band), now stacked per
    TE_overlap_pct bin - so a reader can see whether the "edge vs. internal"
    split holds the same way across the overlap-% range, or whether e.g.
    edge-concentration is mostly a low-overlap-% phenomenon. Stops at 80%:
    classify_overlap_position never computes a 5'/3'/internal/mixed
    breakdown for >=80% ("spans_gene") genes - there's little left to be
    "off to one side of" once the overlap covers most of the transcript -
    so there's no real split to show past that point."""
    cats = list(OVERLAP_POSITION_COLORS)
    df = pc_df[pc_df["TE_overlap_position"].isin(cats)].copy()

    def bin_label(pct):
        for lo, hi in OVERLAP_PCT_BINS:
            if lo <= pct < hi:
                return f"{lo}-{hi}%"
        return None

    df["pct_bin"] = df["TE_overlap_pct"].map(bin_label)
    bin_labels = [f"{lo}-{hi}%" for lo, hi in OVERLAP_PCT_BINS]
    pivot = (
        df.groupby(["pct_bin", "TE_overlap_position"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=bin_labels, columns=cats, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(11, 5.5))
    bottom = pd.Series([0] * len(pivot), index=pivot.index)
    for cat in cats:
        ax.bar(
            pivot.index,
            pivot[cat],
            bottom=bottom,
            label=OVERLAP_POSITION_LABELS[cat],
            color=OVERLAP_POSITION_COLORS[cat],
            zorder=3,
            edgecolor=SURFACE,
            linewidth=1,
        )
        bottom = bottom + pivot[cat]
    for x, total in enumerate(bottom.values):
        ax.annotate(
            f"{int(total)}",
            xy=(x, total),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            color=INK_SECONDARY,
        )
    style_axes(ax)
    ax.set_ylabel("Number of genes")
    ax.set_xlabel("TE overlap (% of transcript) - any 0-80% overlap, not just the 20-80% band")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False, fontsize=8, labelcolor=INK_SECONDARY)
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


FRAGMENT_SIZE_BINS = [(0, 80), (80, 300), (300, 1000), (1000, 3000), (3000, float("inf"))]
FRAGMENT_SIZE_LABELS = ["<80 bp", "80-300 bp", "300-1,000 bp", "1,000-3,000 bp", ">=3,000 bp"]
FRAGMENT_SIZE_COLORS = [CATEGORICAL[7], CATEGORICAL[1], CATEGORICAL[3], CATEGORICAL[0], CATEGORICAL[2]]


FRAGMENT_PCT_BINS = [(0, 20), (20, 80), (80, 100.0001)]


def plot_fragment_size_by_pct(pc_df: pd.DataFrame, title: str, dest: Path):
    """Stacked bar chart, coarser x-axis than plot_overlap_position_by_pct
    (FRAGMENT_PCT_BINS: 0-20/20-50/50-80/80-100, includes the >=80% "likely
    a TE" genes as their own bin), stacked by the size of the largest TE
    element behind the overlap (TE_fragment_max_bp) instead of by overlap
    position - a plain counts-per-bin view instead of a log-log scatter, to
    answer the same question (is a gene's overlap small because the TE
    element itself is small, or because a large, clearly genuine element
    just barely pokes in) more directly. The 80-100% bin dominates the
    y-axis, but that's the point - a reference for what the other bins'
    composition is trending toward."""
    df = pc_df[pc_df["TE_overlap_bp"] > 0].copy()

    def pct_bin_label(pct):
        for lo, hi in FRAGMENT_PCT_BINS:
            if lo <= pct < hi:
                return f"{lo}-{min(hi, 100):g}%"
        return None

    def size_bin_label(bp):
        for (lo, hi), label in zip(FRAGMENT_SIZE_BINS, FRAGMENT_SIZE_LABELS):
            if lo <= bp < hi:
                return label
        return FRAGMENT_SIZE_LABELS[-1]

    df["pct_bin"] = df["TE_overlap_pct"].map(pct_bin_label)
    df["size_bin"] = df["TE_fragment_max_bp"].map(size_bin_label)
    bin_labels = [f"{lo}-{min(hi, 100):g}%" for lo, hi in FRAGMENT_PCT_BINS]
    pivot = (
        df.groupby(["pct_bin", "size_bin"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=bin_labels, columns=FRAGMENT_SIZE_LABELS, fill_value=0)
    )

    fig, ax = plt.subplots(figsize=(11.5, 5.5))
    bottom = pd.Series([0] * len(pivot), index=pivot.index)
    for label, color in zip(FRAGMENT_SIZE_LABELS, FRAGMENT_SIZE_COLORS):
        ax.bar(
            pivot.index, pivot[label], bottom=bottom, label=label, color=color, zorder=3, edgecolor=SURFACE, linewidth=1
        )
        bottom = bottom + pivot[label]
    for x, total in enumerate(bottom.values):
        ax.annotate(
            f"{int(total)}",
            xy=(x, total),
            xytext=(0, 4),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            color=INK_SECONDARY,
        )
    style_axes(ax)
    ax.set_ylabel("Number of genes")
    ax.set_xlabel("TE overlap (% of transcript)")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.legend(
        title="Largest overlapping TE element",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        fontsize=8,
        title_fontsize=8,
        labelcolor=INK_SECONDARY,
    )
    fig.tight_layout()
    fig.savefig(dest, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--te-indir",
        required=True,
        help="Directory containing R_<genome_id>/ RepeatMasker subdirectories "
        "(same layout as summarize_te_content.py --indir)",
    )
    ap.add_argument("--gff-indir", required=True, help="Directory containing one GFF3 file per genome (*.gff / *.gff3)")
    ap.add_argument(
        "--genome-map",
        default=None,
        help="Optional 2-column TSV (genome_id<TAB>gff_path) to override "
        "automatic TE-genome-to-GFF matching, for naming conventions the "
        "automatic prefix-matching can't resolve on its own",
    )
    ap.add_argument("--outdir", required=True, help="Directory to write per-gene TSVs, threshold summaries, and plots")
    ap.add_argument(
        "--min-overlap-bp",
        type=int,
        default=MIN_OVERLAP_BP,
        help=f"A raw TE overlap below this many bp is treated as no overlap at "
        f"all - filters out coincidental boundary-touching noise (default: "
        f"{MIN_OVERLAP_BP})",
    )
    args = ap.parse_args()

    te_indir = Path(args.te_indir)
    gff_indir = Path(args.gff_indir)
    outdir = Path(args.outdir)
    plots_dir = outdir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    genomes = find_genome_dirs(te_indir)
    if not genomes:
        raise SystemExit(f"ERROR: no R_<genome_id>/ directories with .out files found in {te_indir}")

    gff_by_stem = index_gff_dir(gff_indir)
    if not gff_by_stem:
        raise SystemExit(f"ERROR: no .gff/.gff3 files found in {gff_indir}")
    genome_map = load_genome_map(args.genome_map) if args.genome_map else {}

    all_rows = []
    all_transcript_rows = []
    all_skipped_biotypes = Counter()
    gff_seqid_cache: dict = {}
    for genome_id, paths in sorted(genomes.items()):
        te_contigs = {h[0] for h in parse_out(paths["out"])}
        gff_path = find_gff_for_genome(
            genome_id, gff_by_stem, genome_map, te_contigs=te_contigs, seqid_cache=gff_seqid_cache
        )
        if gff_path is None:
            print(
                f"WARNING: no GFF match for TE genome_id {genome_id!r} in "
                f"{gff_indir} ({len(gff_by_stem)} GFF files available) - skipping. "
                f"Use --genome-map to specify this match explicitly."
            )
            continue
        print(f"Processing {genome_id} -> {gff_path.name} ...")

        te_all, te_by_class, te_by_name, element_spans = build_te_tracks(paths["out"])
        gene_meta, rna_to_gene, rna_description, cds_by_rna, exon_by_rna = parse_gff(gff_path)
        rows, transcript_rows, skipped_biotypes = compute_gene_rows(
            genome_id,
            gff_path.stem,
            gene_meta,
            rna_to_gene,
            rna_description,
            cds_by_rna,
            exon_by_rna,
            te_all,
            te_by_class,
            te_by_name,
            element_spans,
            args.min_overlap_bp,
        )
        all_rows.extend(rows)
        all_transcript_rows.extend(transcript_rows)
        all_skipped_biotypes.update(skipped_biotypes)

    if all_skipped_biotypes:
        print(
            f"\nWARNING: {sum(all_skipped_biotypes.values())} gene(s) had a gene_biotype "
            f"that isn't 'protein_coding' or 'pseudogene' (tRNA/rRNA genes included - TE "
            f"overlap is not computed for those) and were skipped entirely (not in "
            f"gene_te_overlap.tsv at all): {dict(all_skipped_biotypes)}"
        )

    # gene_df (every gene, overlap or not) is what n_genes / the threshold
    # tables below are computed from - that denominator has to include the
    # zero-overlap genes to mean anything. The TSV on disk only needs the
    # ones an overlap was actually found for; a ~10,000-row-per-genome file
    # that's 99.7% zeros is not something anyone wants to open.
    gene_df = pd.DataFrame(all_rows)
    nonzero_gene_df = gene_df[gene_df["TE_overlap_bp"] > 0]
    # TE_fragment_max_bp* describe the single largest-by-its-own-size TE
    # element touching the gene, not the gene's total overlap (that's
    # TE_overlap_bp/TE_overlap_pct) - reading them side by side in the TSV
    # reads as contradictory (a "29 bp" max-fragment-overlap next to an
    # 80%+ TE_overlap_pct looks like an error, when it just means the
    # biggest single element only clips the gene and other fragments make
    # up the rest). Kept internally for plot_fragment_size_by_pct below,
    # dropped from the file itself.
    _FRAGMENT_MAX_COLS = ["TE_fragment_max_bp", "TE_fragment_max_bp_name", "TE_fragment_max_bp_overlap_bp"]
    nonzero_gene_df.drop(columns=_FRAGMENT_MAX_COLS).to_csv(outdir / "gene_te_overlap.tsv", sep="\t", index=False)
    print(
        f"\nWrote {outdir / 'gene_te_overlap.tsv'} "
        f"({len(nonzero_gene_df)} genes with TE overlap, out of {len(gene_df)} total)"
    )

    transcript_df = pd.DataFrame(all_transcript_rows)
    transcript_df.to_csv(outdir / "multi_transcript_te_overlap.tsv", sep="\t", index=False)
    n_multi_overlap_genes = transcript_df["gene_id"].nunique() if len(transcript_df) else 0
    print(
        f"Wrote {outdir / 'multi_transcript_te_overlap.tsv'} "
        f"({n_multi_overlap_genes} multi-transcript genes with TE overlap, "
        f"{len(transcript_df)} transcript rows)"
    )

    pc_df = gene_df[gene_df["gene_category"] == "protein_coding"]
    pg_df = gene_df[gene_df["gene_category"] == "pseudogene"]
    print(f"  protein_coding: {len(pc_df)}   pseudogene: {len(pg_df)}")

    pc_thresh = threshold_table(pc_df)
    pc_thresh.to_csv(outdir / "protein_coding_threshold_summary.tsv", sep="\t", index=False)
    pg_thresh = threshold_table(pg_df)
    pg_thresh.to_csv(outdir / "pseudogene_threshold_summary.tsv", sep="\t", index=False)
    print(f"Wrote {outdir / 'protein_coding_threshold_summary.tsv'}")
    print(f"Wrote {outdir / 'pseudogene_threshold_summary.tsv'}")

    # --- Plots ---------------------------------------------------------
    plot_hist(
        pc_df["TE_overlap_pct"],
        "Protein-coding genes: CDS overlap with TEs (all genes)",
        plots_dir / "protein_coding_te_overlap_hist.png",
        log_y=True,
    )
    pc_nonzero_df = pc_df[pc_df["TE_overlap_bp"] > 0]
    plot_hist(
        pc_nonzero_df["TE_overlap_pct"],
        "Protein-coding genes with TE overlap > 0",
        plots_dir / "protein_coding_te_overlap_hist_nonzero.png",
        log_y=False,
        color=CATEGORICAL[1],
    )
    plot_hist_by_genome(
        pc_nonzero_df,
        "Protein-coding CDS/TE overlap by genome (genes with overlap > 0)",
        plots_dir / "protein_coding_te_overlap_by_genome.png",
        bins=20,
        log_y=False,
    )
    pc_total_row = pc_thresh[pc_thresh["genome_id"] == "ALL"].iloc[0]
    plot_threshold_heatmap(
        pc_thresh, "TE-overlap thresholds by genome", plots_dir / "protein_coding_threshold_sensitivity_by_genome.png"
    )
    plot_band_bar(
        pc_thresh,
        "20-80%",
        "Genes in the 20-80% band",
        "Genes in the 20-80% 'inspect' band, by genome",
        plots_dir / "protein_coding_inspect_band_by_genome.png",
    )
    print(
        f"\n{int(pc_total_row['20-80%'])} protein-coding genes fall in the 20-80% "
        f"'inspect' band (pooled across genomes) - see the '20-80%' column in "
        f"protein_coding_threshold_summary.tsv for the per-genome counts."
    )

    plot_band_bar(
        pc_thresh,
        ">=80%",
        "Genes in the >=80% band",
        "Genes in the >=80% 'eliminated' band, by genome",
        plots_dir / "protein_coding_eliminated_band_by_genome.png",
        color=CATEGORICAL[2],
    )
    print(
        f"{int(pc_total_row['>=80%'])} protein-coding genes fall in the >=80% "
        f"'eliminated' band (pooled across genomes) - see the '>=80%' column in "
        f"protein_coding_threshold_summary.tsv for the per-genome counts."
    )

    plot_overlap_position_bar(
        pc_df,
        "Where partial TE overlaps fall along the gene (20-80% band)",
        plots_dir / "protein_coding_overlap_position.png",
    )
    plot_overlap_position_by_pct(
        pc_df,
        "Overlap position by TE-overlap % (0-80% overlap)",
        plots_dir / "protein_coding_overlap_position_by_pct.png",
    )
    plot_fragment_size_by_pct(
        pc_df,
        "Size of the TE element behind the overlap, by TE-overlap %",
        plots_dir / "protein_coding_te_fragment_size.png",
    )
    n_check_model = int(pc_df["check_gene_model"].sum())
    print(
        f"{n_check_model} protein-coding genes have TE overlap concentrated at one end "
        f"(5' or 3') of the gene - candidates for a gene-model check, not just a TE-content "
        f"call. See the 'TE_overlap_position'/'check_gene_model' columns in gene_te_overlap.tsv."
    )

    plot_hist(
        pg_df["TE_overlap_pct"],
        "Pseudogene overlap with TEs (pooled)",
        plots_dir / "pseudogene_te_overlap_hist.png",
        log_y=True,
        color=CATEGORICAL[2],
    )

    if len(transcript_df):
        per_gene = transcript_df.groupby(["genome_id", "gene_id"])["transcript_TE_overlap_pct"].agg(["min", "max"])
        # "divergent" = isoforms disagree a lot (one mostly TE-free, another
        # mostly TE) - a >=20 percentage-point gap between a gene's lowest-
        # and highest-overlap isoform. A gene where every isoform sits near
        # the same value (all ~100%, or all ~5%) is NOT divergent, even if
        # its min happens to be low - that just means the gene as a whole
        # has low overlap, not that isoform choice matters.
        divergent = per_gene[(per_gene["max"] - per_gene["min"]) >= 20]
        print(
            f"\n{len(per_gene)} multi-transcript protein_coding genes have TE overlap in at "
            f"least one isoform; {len(divergent)} of those show a real split between isoforms "
            f"(>=20 percentage points between the lowest- and highest-overlap isoform) - see "
            f"multi_transcript_te_overlap.tsv."
        )

    print(f"\nPlots written to {plots_dir}/")


if __name__ == "__main__":
    main()
