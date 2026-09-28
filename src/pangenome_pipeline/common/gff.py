"""Finding and reading GFF files."""

import os
from pathlib import Path
from urllib.parse import unquote_plus

# Gene-level feature types (GFF3 column 3) mapped to the biotype they imply when the record
# carries no gene_biotype/biotype/gene_type attribute at all. NCBI's annotwriter always sets
# gene_biotype=... explicitly on a plain "gene" feature, so that attribute wins when present
# (see infer_biotype). Sources that instead encode biotype in the feature type itself -
# VEuPathDB (FungiDB, PlasmoDB, ...) and Ensembl-style GFF3 both do this - are the ones this
# table is for: a bare "gene" with no biotype attribute is, by GFF3 convention, protein-coding;
# RNA genes get their own explicit type.
GENE_FTYPE_TO_BIOTYPE = {
    "gene": "protein_coding",
    "protein_coding_gene": "protein_coding",
    "pseudogene": "pseudogene",
    "tRNA_gene": "tRNA",
    "rRNA_gene": "rRNA",
    "ncRNA_gene": "ncRNA",
    "snoRNA_gene": "snoRNA",
    "snRNA_gene": "snRNA",
    "miRNA_gene": "miRNA",
}
BIOTYPE_ATTR_KEYS = ("gene_biotype", "biotype", "gene_type")

# gene_biotype=/biotype= attribute values (NCBI/GenBank and Ensembl-style conventions),
# lower-cased, to the broad category the pipeline reports and filters on.
BIOTYPE_TO_CATEGORY = {
    "protein_coding": "protein_coding",
    "mrna": "protein_coding",
    "pseudogene": "pseudogene",
    "trna": "tRNA",
    "rrna": "rRNA",
    "ncrna": "ncRNA",
    "lncrna": "ncRNA",
    "lnc_rna": "ncRNA",
    "snrna": "ncRNA",
    "snorna": "ncRNA",
    "mirna": "ncRNA",
    "misc_rna": "ncRNA",
}

# Free-text description, checked in this order - NCBI puts "product" on the mRNA/CDS rather
# than the gene; VEuPathDB (FungiDB, ...) puts "description" directly on the gene; "Note" and
# "Name" are the common fallbacks when neither of those is set.
DESCRIPTION_ATTR_KEYS = ("description", "product", "Note", "note", "Name")


def find_gff(gff_dir, name, exts=(".gff3", ".gff"), case_insensitive=False):
    """Path of <gff_dir>/<name><ext> for the first ext that exists, or None.

    With case_insensitive, falls back to any *.gff* file whose stem matches name ignoring case.
    Returns a Path if gff_dir is a Path, otherwise a string, so callers keep the type they use.
    """
    as_path = isinstance(gff_dir, Path)
    for ext in exts:
        candidate = os.path.join(gff_dir, f"{name}{ext}")
        if os.path.isfile(candidate):
            return Path(candidate) if as_path else candidate
    if case_insensitive:
        lowered = name.lower()
        for f in Path(gff_dir).glob("*.gff*"):
            if f.stem.lower() == lowered:
                return f if as_path else str(f)
    return None


def gff_seqids(path):
    """Every column-1 sequence name in a GFF (comment and blank lines skipped)."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ids.add(line.split("\t", 1)[0])
    return ids


def infer_biotype(ftype, attrs):
    """The record's own gene_biotype/biotype/gene_type attribute if it has one (NCBI/RefSeq,
    Ensembl-style, GENCODE - that explicit label is authoritative), otherwise what the feature
    type implies (GENE_FTYPE_TO_BIOTYPE), otherwise "NA"."""
    for key in BIOTYPE_ATTR_KEYS:
        if key in attrs:
            return attrs[key]
    return GENE_FTYPE_TO_BIOTYPE.get(ftype, "NA")


def biotype_category(value, other_rna_label="ncRNA"):
    """Category for a gene_biotype=/biotype= value: exact matches first, then the same idea for
    the many variants other sources use (processed_pseudogene, Mt_tRNA, rRNA_pseudogene, sRNA,
    scaRNA, lncRNA...). Anything else (e.g. IG_C_gene) is None - left for the caller to report.

    RNA genes other than tRNA/rRNA come back as other_rna_label ("ncRNA" by default)."""
    v = value.strip().strip('"').lower()
    category = BIOTYPE_TO_CATEGORY.get(v)
    if category is None:
        if "pseudogene" in v:
            category = "pseudogene"
        elif "trna" in v:
            category = "tRNA"
        elif "rrna" in v:
            category = "rRNA"
        elif "rna" in v:
            category = "ncRNA"
    return other_rna_label if category == "ncRNA" else category


def parse_attributes(attr_str):
    """GFF3 column 9 as a dict; values are left URL-encoded."""
    attrs = {}
    for kv in attr_str.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            attrs[k] = v
    return attrs


def gene_description(attrs):
    """The first non-empty description-like attribute, decoded: %XX escapes and '+' as a space
    (VEuPathDB writes "hypothetical+protein"). "" if there is none."""
    for key in DESCRIPTION_ATTR_KEYS:
        v = attrs.get(key)
        if v:
            return unquote_plus(v)
    return ""
