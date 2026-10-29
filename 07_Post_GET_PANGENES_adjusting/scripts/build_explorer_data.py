#!/usr/bin/env python3
"""
Build the data files for the interactive pangenome explorer (explorer.html), from a get_pangenes matrix and the
matching GFFs (build_explorer_job.sh runs it, on the get_pangenes matrix after step 06 or on the final pangenome
after step 07). Gene quality comes from step 05's orf_status: High (valid, repaired_miniprot, selenoprotein), Low
(invalid, pseudogene) or Not judged (ORF broken by the assembly: invalid_assembly_gap, invalid_near_contig_end); the
pangene quality is the majority of its judged genes (High / Low / Mixed / Not judged).

With --absences (audit_absences.py's absences.tsv, or step 07's absences_explained.tsv, for this same pangenome), every pangene also gets "absent": the
category of every genome where it has no gene (1a ... 8), and "unknown": the genomes where the absence comes from the
assembly (categories 6-8), for pangenes in >= --flag-min-genomes genomes; the explorer shows these as "?".
--absence-flags (step 07's absence_flags.tsv) gives only "unknown".

With --winners (step 05's winners.tsv), every pangene also gets "panel_genes": the step 05 panel genes (basket=) of its
genes, each with the PSAURON score of its winner model and how many genes of the pangene it has, shown once per pangene.

For every pangene, for every gene it contains, collects step 05's attributes (ID, original_ID, model_from,
model_origin, basket, orf_status, description, ebi_biotype, has_domain, tandem_array_copy,
not_in_reference_annotation, pseudogene, coverage, sequence_ID - whichever the gene has), coordinates, and its
InterPro/Pfam domain ids (from CDS Dbxref, name looked up in --entry-list when given). The "source" field is
model_origin (where the model came from).

Output is one JSON object, kept small with a columnar layout (a "genes" field list + one array
per gene, not one object per gene - about half the size of the naive encoding) rather than by
lowering how much detail is in it:
  {
    "genomes": [genome names, in matrix column order],
    "gene_fields": [attribute names, fixed order],
    "pangenes": [
      {"id": pangene_id, "occupancy": N, "category": "Core"|"Soft-core"|"Shell"|"Cloud",
       "genes": {genome_name: [[gene_field_values...], ...one array per copy in that genome]}}
      , ...
    ]
  }

Usage:
    python3 build_explorer_data.py \
        --genes-matrix pangenome_final/pangene_matrix_genes.tr.tab --gff-dir pangenome_final/gff \
        --absences pangenome_final/absences_explained.tsv \
        --entry-list entry.list -o pangenome_explorer/data \
        [--core-missing 0] [--softcore-missing 1] [--cloud-max 2]

Needs only the standard library.
"""
import argparse
import csv
import glob
import json
import os
import re
import sys
from collections import Counter, defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
DBXREF_RE = re.compile(r"(?:^|;)Dbxref=([^;]+)")

GENE_ATTRS = [
    "ID", "original_ID", "model_from", "model_origin", "basket", "orf_status", "description", "ebi_biotype",
    "has_domain", "tandem_array_copy", "not_in_reference_annotation", "pseudogene", "coverage", "sequence_ID",
]
# fixed record layout written for every gene copy, in this order
GENE_FIELDS = ["source"] + GENE_ATTRS + ["contig", "start", "end", "strand", "interpro", "pfam", "quality_class"]

# ---- ORF-based gene quality, from step 05's orf_status ----
HIGH = {"valid", "repaired_miniprot", "selenoprotein"}
LOW = {"invalid", "pseudogene"}


def classify_gene(orf_status):
    s = (orf_status or "").strip()
    return "High" if s in HIGH else "Low" if s in LOW else "Not judged"


def pangene_quality(classes):
    judged = [c for c in classes if c != "Not judged"]
    if not judged:
        return "Not judged"
    n_high = sum(1 for c in judged if c == "High")
    n_low = len(judged) - n_high
    if n_high > n_low:
        return "High"
    if n_low > n_high:
        return "Low"
    return "Mixed"


def parse_attrs(field9):
    out = {}
    for kv in field9.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def categorize(occ, n, core_missing=0, softcore_missing=1, cloud_max=2):
    if occ >= n - core_missing:
        return "Core"
    if occ >= n - softcore_missing:
        return "Soft-core"
    if occ <= cloud_max:
        return "Cloud"
    return "Shell"


def load_contig_sizes(fai_dir, genome_names, suffix=".fai"):
    """genome -> {contig: length}, from samtools-faidx .fai files (col 1 = name, col 2 = length)."""
    out = {}
    for g in genome_names:
        path = os.path.join(fai_dir, f"{g}{suffix}")
        sizes = {}
        if os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 2:
                        sizes[f[0]] = int(f[1])
        else:
            print(f"WARNING: no {path}, contig lengths unavailable for {g} (quality classification treats its genes as not-short)", file=sys.stderr)
        out[g] = sizes
    return out


def load_entry_names(path):
    names = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[0].startswith("IPR"):
                    names[f[0]] = f[2]
    return names


def scan_gff(path):
    """gene_id -> record dict (attrs + coords + source), gene_id -> {interpro:set, pfam:set}."""
    genes = {}
    mrna_gene = {}
    domains = defaultdict(lambda: {"interpro": set(), "pfam": set()})
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            contig, src, ftype, start, end, strand, attrs = f[0], f[1], f[2], f[3], f[4], f[6], f[8]
            if ftype == "gene" or ftype == "pseudogene":
                a = parse_attrs(attrs)
                gid = a.get("ID")
                if not gid:
                    continue
                genes[gid] = {"source": src, "contig": contig, "start": start, "end": end,
                               "strand": strand, "attrs": a}
            elif ftype in ("mRNA", "transcript"):
                a = parse_attrs(attrs)
                i, p = a.get("ID"), a.get("Parent")
                if i and p:
                    mrna_gene[i] = p
            elif ftype == "CDS":
                p, d = PARENT_RE.search(attrs), DBXREF_RE.search(attrs)
                if p and d:
                    gid = mrna_gene.get(p.group(1), p.group(1))
                    for t in d.group(1).split(","):
                        if t.startswith("InterPro:"):
                            domains[gid]["interpro"].add(t.split(":", 1)[1])
                        elif t.startswith("Pfam:"):
                            domains[gid]["pfam"].add(t.split(":", 1)[1])
    return genes, domains


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-matrix", required=True)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--gff-suffix", default="_basket_annotation.gff")
    ap.add_argument("--entry-list", default="")
    ap.add_argument("--contig-sizes-dir", default="", help="directory with <genome>.fai (samtools faidx); used for the provisional quality classification's short-contig check")
    ap.add_argument("--core-missing", type=int, default=0)
    ap.add_argument("--softcore-missing", type=int, default=1)
    ap.add_argument("--cloud-max", type=int, default=2)
    ap.add_argument("--absences", default="", help="audit_absences.py's absences.tsv for this pangenome: category of every absence")
    ap.add_argument("--flag-min-genomes", type=int, default=3, help="with --absences: '?' only for pangenes in at least this many genomes")
    ap.add_argument("--winners", default="", help="step 05's winners.tsv: PSAURON of every panel gene (shown per pangene)")
    ap.add_argument("--absence-flags", default="", help="step 07's absence_flags.tsv (pangene, genome, code): shown as '?'")
    ap.add_argument("-o", "--out", required=True, help="output directory for manifest.json + chunk_*.json")
    ap.add_argument("--chunk-mb", type=float, default=12.0, help="target size per chunk file, MB (default 12, stays under the 16MB per-file publish limit)")
    args = ap.parse_args()

    entry_name = load_entry_names(args.entry_list)

    with open(args.genes_matrix, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header = rows[0]
    genome_cols = [i for i, g in enumerate(header) if i > 0 and g.strip()]
    genome_names = [header[i].split("_for_get_pangenes")[0] for i in genome_cols]
    n = len(genome_cols)

    print(f"Scanning {n} genome GFF(s)...", file=sys.stderr)
    gff_genes, gff_domains = {}, {}
    for gname in genome_names:
        path = os.path.join(args.gff_dir, f"{gname}{args.gff_suffix}")
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")
        genes, domains = scan_gff(path)
        gff_genes[gname] = genes
        gff_domains[gname] = domains
        print(f"  {gname}: {len(genes)} genes", file=sys.stderr)

    contig_sizes = load_contig_sizes(args.contig_sizes_dir, genome_names) if args.contig_sizes_dir else {g: {} for g in genome_names}

    def gene_record(gname, gid):
        g = gff_genes[gname].get(gid)
        if g is None:
            return None
        a = g["attrs"]
        d = gff_domains[gname].get(gid, {"interpro": set(), "pfam": set()})
        ipr = sorted(d["interpro"])
        ipr_named = [f"{i} {entry_name[i]}" if i in entry_name else i for i in ipr]
        has_domain = bool(ipr) or bool(d["pfam"])
        contig_len = contig_sizes.get(gname, {}).get(g["contig"])
        qclass = classify_gene(a.get("orf_status"))
        rec = [a.get("model_origin") or g["source"]]
        for k in GENE_ATTRS:
            rec.append(a.get(k, ""))
        rec += [g["contig"], g["start"], g["end"], g["strand"], ipr_named, sorted(d["pfam"]), qclass]
        return rec

    unknown, absent_cat, absent_place = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    if args.absences:
        with open(args.absences, newline="") as fh:
            for fr in csv.DictReader(fh, delimiter="\t"):
                absent_cat[fr["pangene"]][fr["genome"]] = fr.get("code") or fr.get("audit_code")
                if fr.get("place"):
                    absent_place[fr["pangene"]][fr["genome"]] = fr["place"]   # contig:start-end, for the locus view
    if args.absence_flags:
        with open(args.absence_flags, newline="") as fh:
            for fr in csv.DictReader(fh, delimiter="\t"):
                unknown[fr["pangene"]][fr["genome"]] = fr["code"]

    psauron = {}
    if args.winners:
        with open(args.winners, newline="") as fh:
            for fr in csv.DictReader(fh, delimiter="\t"):
                psauron[fr["gene"]] = fr.get("psauron", "")
    bi = GENE_FIELDS.index("basket")

    pangenes_out = []
    n_missing = 0
    for r in rows[1:]:
        if not r or r[0].startswith("chr"):
            continue
        pangene_id = r[0]
        genes_by_genome = {}
        occ = 0
        for gi, gname in zip(genome_cols, genome_names):
            if gi >= len(r) or r[gi] in ("", "-", "NA"):
                continue
            gene_ids = r[gi].split(",")
            recs = []
            for gid in gene_ids:
                rec = gene_record(gname, gid)
                if rec is None:
                    n_missing += 1
                    continue
                recs.append(rec)
            if recs:
                genes_by_genome[gname] = recs
                occ += 1
        cat = categorize(occ, n, args.core_missing, args.softcore_missing, args.cloud_max)
        qmajority = pangene_quality([rec[-1] for recs in genes_by_genome.values() for rec in recs])  # quality_class is the last field
        rec_p = {"id": pangene_id, "occupancy": occ, "category": cat,
                 "quality_majority": qmajority, "genes": genes_by_genome}
        if args.winners:
            nb = Counter(rec[bi] for recs in genes_by_genome.values() for rec in recs if rec[bi])
            rec_p["panel_genes"] = [[b, psauron.get(b, ""), k] for b, k in nb.most_common()]
        ab = {g: c for g, c in absent_cat.get(pangene_id, {}).items() if g not in genes_by_genome}
        if ab:
            rec_p["absent"] = ab
            pl = {g: absent_place[pangene_id][g] for g in ab if g in absent_place.get(pangene_id, {})}
            if pl:
                rec_p["absent_place"] = pl
        unk = {g: c for g, c in unknown.get(pangene_id, {}).items() if g not in genes_by_genome}
        if not unk and occ >= args.flag_min_genomes:
            unk = {g: c for g, c in ab.items() if c in ("6", "7", "8")}
        if unk:
            rec_p["unknown"] = unk
        pangenes_out.append(rec_p)

    if n_missing:
        print(f"WARNING: {n_missing} gene id(s) from the matrix were not found in their genome's GFF", file=sys.stderr)

    os.makedirs(args.out, exist_ok=True)

    # Split pangenes_out into chunk files of roughly --chunk-mb each, so every published file
    # stays comfortably under the 16MB per-file publish limit. Sized from a quick sample rather
    # than by pangene count, since gene density per pangene varies a lot (Core vs Cloud).
    sample_n = min(200, len(pangenes_out))
    sample_bytes = len(json.dumps(pangenes_out[:sample_n], separators=(",", ":")))
    bytes_per_pangene = sample_bytes / max(sample_n, 1)
    per_chunk = max(1, int((args.chunk_mb * 1024 * 1024) / bytes_per_pangene))

    chunk_files = []
    for start in range(0, len(pangenes_out), per_chunk):
        chunk = pangenes_out[start:start + per_chunk]
        fname = f"chunk_{len(chunk_files):03d}.json"
        with open(os.path.join(args.out, fname), "w") as fh:
            json.dump(chunk, fh, separators=(",", ":"))
        chunk_files.append(fname)

    manifest = {"genomes": genome_names, "gene_fields": GENE_FIELDS, "chunks": chunk_files,
                "n_pangenes": len(pangenes_out)}
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, separators=(",", ":"))

    total_mb = sum(os.path.getsize(os.path.join(args.out, f)) for f in chunk_files) / (1024 * 1024)
    print(f"{len(pangenes_out)} pangenes written across {len(chunk_files)} chunk(s) "
          f"({total_mb:.2f} MB total) + manifest.json, in {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
