"""Parsing and overlap logic shared by the step 05 scripts that read <genome>_final_annotation.gff:
filter_final_annotation.py (7a), tag_final_annotation.py (7b) and the two QC plot scripts.

Not run on its own - imported by those scripts, which sit in the same folder.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote_plus

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from pangenome_common.intervals import segments_overlap_bp

ID_RE = re.compile(r"ID=([^;]+)")
PARENT_RE = re.compile(r"Parent=([^;]+)")
DESCRIPTION_RE = re.compile(r"(?:^|;)description=([^;]+)")
PRODUCT_RE = re.compile(r"(?:^|;)product=([^;]+)")

# Feature types that only ever appear nested under a level-1 gene feature, never as an
# independent top-level entity - used to recognize "a new gene-level block just started"
# in the output-writing passes, without having to enumerate every naming convention
# an annotation might use for its non-protein-coding gene types (NCBI's "gene", FungiDB's
# "ncRNA_gene", etc.) - anything not in this set is treated as a new top-level block.
CHILD_FEATURE_TYPES = {
    "mRNA",
    "CDS",
    "exon",
    "five_prime_UTR",
    "three_prime_UTR",
    "ncRNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "pseudogenic_transcript",
}


def int_attr(attrs, key):
    m = re.search(rf"{key}=(\d+)", attrs)
    return int(m.group(1)) if m else None


def float_attr(attrs, key):
    m = re.search(rf"{key}=([\d.]+)", attrs)
    return float(m.group(1)) if m else None


def valid_orf(attrs):
    """True/False from Liftoff's valid_ORF= (or valid_ORFs=N >= 1); None if neither is present."""
    if "valid_ORF=True" in attrs:
        return True
    if "valid_ORF=False" in attrs:
        return False
    m = re.search(r"valid_ORFs=(\d+)", attrs)
    if m:
        return int(m.group(1)) >= 1
    return None


def parse_gff(path):
    """Gene-level records, mRNA->gene links, and per-mRNA CDS segments/domain flags."""
    genes = {}
    mrna_to_gene = {}
    cds_by_mrna = defaultdict(list)
    mrna_has_domain = defaultdict(bool)

    with open(path) as fh:
        raw_lines = fh.readlines()

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        seqname, source, feature, start, end, score, strand = fields[:7]
        attrs = fields[8]
        m = ID_RE.search(attrs)
        feat_id = m.group(1) if m else None

        if feature == "protein_coding_gene" and feat_id:
            desc_m = DESCRIPTION_RE.search(attrs) or PRODUCT_RE.search(attrs)
            description = unquote_plus(desc_m.group(1).strip().lower()) if desc_m else ""
            genes[feat_id] = {
                "id": feat_id,
                "seqname": seqname,
                "source": source,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "extra_copy_number": int_attr(attrs, "extra_copy_number"),
                "valid_ORF": valid_orf(attrs),
                "coverage": float_attr(attrs, "coverage"),
                "description": description,
            }
        elif feature == "mRNA" and feat_id:
            parent = PARENT_RE.search(attrs)
            if parent:
                mrna_to_gene[feat_id] = parent.group(1)

    for line in raw_lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        feature, start, end, attrs = fields[2], fields[3], fields[4], fields[8]
        if feature != "CDS":
            continue
        parent = PARENT_RE.search(attrs)
        if not parent:
            continue
        mrna_id = parent.group(1)
        if mrna_to_gene.get(mrna_id) not in genes:
            continue
        cds_by_mrna[mrna_id].append((int(start), int(end)))
        if "InterPro:" in attrs or "PFAM:" in attrs or "Pfam:" in attrs:
            mrna_has_domain[mrna_id] = True

    for gene_id, gene in genes.items():
        gene_mrnas = [m for m, g in mrna_to_gene.items() if g == gene_id]
        best_mrna, best_bp = None, -1
        for mrna_id in gene_mrnas:
            segs = cds_by_mrna.get(mrna_id, [])
            bp = sum(e - s + 1 for s, e in segs)
            if bp > best_bp:
                best_bp, best_mrna = bp, mrna_id
        # sorted left-to-right by genomic position - needed by classify_minor_position below,
        # which walks segments outermost-to-innermost.
        gene["cds_segments"] = sorted(cds_by_mrna.get(best_mrna, []))
        gene["cds_bp"] = best_bp if best_bp > 0 else 0
        gene["has_domain"] = any(mrna_has_domain.get(m, False) for m in gene_mrnas)

    return genes


def find_cds_overlaps(genes, min_overlap_bp=1, with_strand=False):
    """Sweep-line scan for CDS-to-CDS overlapping gene pairs: (gene1_id, gene2_id, overlap_bp),
    plus a same-strand flag as a 4th item with with_strand. Pairs overlapping by less than
    min_overlap_bp are left out (default 1: any real overlap counts)."""
    by_seq = defaultdict(list)
    for gene in genes.values():
        by_seq[gene["seqname"]].append(gene)

    pairs = []
    for seqname, seq_genes in by_seq.items():
        seq_genes.sort(key=lambda g: g["start"])
        active = []
        for gene in seq_genes:
            active = [a for a in active if a["end"] >= gene["start"]]
            for other in active:
                ov_bp = segments_overlap_bp(other["cds_segments"], gene["cds_segments"])
                if ov_bp >= min_overlap_bp:
                    pair = (other["id"], gene["id"], ov_bp)
                    pairs.append((*pair, other["strand"] == gene["strand"]) if with_strand else pair)
            active.append(gene)
    return pairs


def classify_minor_position(bigger_gene, smaller_cds_segments):
    """Where a minor overlap sits on the bigger gene's CDS: "flank" (only touches its first or
    last CDS segment - the two gene models plausibly disagree only at the boundary, both genes
    likely real), "inside" (touches an internal segment - more likely the smaller gene is a
    spurious call sitting in what's actually an intron), or "single_exon" (bigger gene has only
    one CDS segment total, flank/inside doesn't apply)."""
    segs = bigger_gene["cds_segments"]
    if len(segs) <= 1:
        return "single_exon"
    outer = {segs[0], segs[-1]}
    overlapped = {s for s in segs for ss, se in smaller_cds_segments if max(s[0], ss) <= min(s[1], se)}
    if not overlapped:
        return "inside"  # shouldn't happen (ov_bp>0 implies some segment matched) - fall through safely
    return "flank" if overlapped <= outer else "inside"


def assign_new_ids(genes, gene_ids, strain_index, id_prefix):
    """gene_id -> "<id_prefix><strain_index>g<00001...>" for gene_ids, numbered by position.

    Genes starting at the same position (e.g. Liftoff extra copies) are ordered by their original
    ID, so the numbering is the same on every run whatever order gene_ids comes in."""
    ordered = sorted(gene_ids, key=lambda gid: (genes[gid]["seqname"], genes[gid]["start"], gid))
    width = max(5, len(str(len(ordered))))
    return {gid: f"{id_prefix}{strain_index}g{i:0{width}d}" for i, gid in enumerate(ordered, start=1)}


def add_original_id_attr(attrs, old_id):
    return attrs.rstrip(";") + f";original_ID={old_id}"


def source_counts(genes, reference_strain):
    """reference_strain is the actual genome name (e.g. "AF293"), not the "main" tag - needed
    to split "rescued_<reference>" out from every other strain's "rescued_<x>"."""
    ref_rescued_tag = f"rescued_{reference_strain}"
    counts = {"main": 0, "original": 0, "reference_rescued": 0, "good_others": 0, "rescued_others": 0}
    for g in genes.values():
        src = g["source"]
        if src == "main":
            counts["main"] += 1
        elif src == ref_rescued_tag:
            counts["reference_rescued"] += 1
        elif src.startswith("good_"):
            counts["good_others"] += 1
        elif src.startswith("rescued_"):
            counts["rescued_others"] += 1
        else:
            counts["original"] += 1
    return counts
