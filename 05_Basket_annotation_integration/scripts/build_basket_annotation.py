#!/usr/bin/env python3
"""
Write each genome's final GFF from the basket decisions (one gene per panel gene and genome), ready for get_pangenes.

Inputs:
  baskets_v2.tsv          the baskets and their winners (baskets_v2.py); one winner = one panel gene (a basket without
                          a winner writes nothing, and its genes do not hold their place)
  overlap votes           optional (--overlap-votes, basket_overlaps.py): genes that lose an overlap, removed from every
                          genome. The basket job no longer uses it: the overlap rule runs at the end (overlap_filter.py)
  decisions.tsv           what each genome takes for each panel gene, in both modes (decide_models.py)
  native GFFs, Liftoff transfers, the reference GFF before its curation (to find the genes the curation removed)

Removed panel genes:
  overlap   (only with --overlap-votes) the gene loses an overlap (pair rules) in at least one genome
  curation  a panel gene without a reference gene whose model, in the reference genome, lies >= 50% over a gene that
            the curation of the reference removed (present in --old-reference-gff, absent from the reference GFF; or,
            with --curation-removed-list, only the genes of that list, e.g. when the curation also renamed genes)
Per genome (--mode keep_similar | prefer_panel_model, column decision_<mode> of decisions.tsv):
  keep_original                   (keep_similar) the genome's own gene is written unchanged (model_origin=original)
  target (own model)              the panel gene in its own genome (model_origin=original)
  replace / fill / target_invalid the panel gene's Liftoff transfer (model_origin=replaced / filled / target_invalid_orf)
  alternative                     another member's model with a valid ORF (model_origin=alternative)
  rescued                         a valid model covering less of the panel gene (model_origin=rescued)
  absent                          nothing
A transfer is used only at the panel gene's own place: when it lies >= 50% over a native gene of another basket, it is
not used there (that place belongs to the other basket) and the genome counts as absent for this panel gene. Genes
of a basket removed as false (overlap, curation) do not hold their place.

Output GFF (one per genome): one gene, one transcript, its exons and CDS. Every gene gets a new unique ID,
<id-prefix><genome index>g<number> (genome index from --genomes, numbered along the genome), and in column 9:
  original_ID   the genome's own gene ID when its gene is kept or replaced (empty when filled)
  model_from    the gene the model comes from
  model_origin  original / replaced / filled / alternative / rescued / target_invalid_orf
  basket        the panel gene (winner) ID
  has_domain    InterPro/Pfam domain of model_from
  not_in_reference_annotation=True   panel genes without a reference gene
  tandem_array_copy=True             a copy of a tandem array (decide_models.py's array_models.tsv)
  replaced_broken_model              the array's model(s) with a broken ORF that a genome's own valid copy replaced
  protected=True                     a gene of a protected basket (the shield, --protected-baskets): it wins
                                     "one place, one gene" against unprotected models of non-reference baskets
  description   kept from the source gene
Genes left without a model (--members, baskets_v2.py's members.tsv): after the steps above, a member of a basket
with a winner that no model of its basket covers is
  - a piece of the gene, when its basket has a model in this genome (copies were split off in baskets_v2.py): left out;
  - otherwise (the basket's transfer was not used, or was dropped by "one place, one gene"): it keeps its own model at
    its own place (model_origin=original_other_place) when it passes members.tsv's quality check and does not lie
    >= 50% (same strand) over a model already written. With several such members, the one with the best PSAURON
    (then the largest protein) comes back; the others are pieces of it.
Every member left out is listed in left_out.tsv with the reason. No member is left out silently.
The models dropped by "one place, one gene" (same strand, >= 50% of the smaller) are written to
<genome>_stacked_dropped.gff (IDs <id-prefix><genome index>s<number>), so that overlap_filter.py can put one back when
the gene that took its place is removed later.
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
DESC_RE = re.compile(r"(?:^|;)description=([^;]+)")


def gene_blocks(path, want, primary_only=False):
    """gene ID -> (gene fields, longest transcript's exon and CDS fields) for the genes in want."""
    genes, tx_of, kids = {}, {}, defaultdict(list)
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        i, p = ID_RE.search(f[8]), PARENT_RE.search(f[8])
        if i and not p:
            if i.group(1) in want and not (primary_only and "extra_copy_number=" in f[8]
                                           and "extra_copy_number=0" not in f[8]):
                genes[i.group(1)] = f
        elif p:
            par = p.group(1).split(",")[0]
            if f[2] in ("mRNA", "transcript"):
                tx_of[i.group(1) if i else ""] = par
            kids[par].append(f)
    out = {}
    for g, gf in genes.items():
        best, best_len = None, -1
        for t in kids.get(g, []):
            if t[2] not in ("mRNA", "transcript") or t[0] != gf[0]:
                continue
            if not (int(t[3]) >= int(gf[3]) and int(t[4]) <= int(gf[4])):
                continue
            ti = ID_RE.search(t[8]).group(1)
            cds = [c for c in kids.get(ti, []) if c[2] == "CDS"]
            L = sum(int(c[4]) - int(c[3]) + 1 for c in cds)
            if L > best_len:
                best, best_len = (t, [c for c in kids.get(ti, []) if c[2] in ("exon", "CDS")]), L
        if best and best_len > 0:
            out[g] = (gf, best[0], best[1])
    return out


def keep_attrs(col9, drop=()):
    """the attributes of a source line, without ID, Parent and the keys in drop (they are set again)."""
    out = []
    for kv in col9.rstrip("\n").split(";"):
        key = kv.split("=", 1)[0]
        if kv and key not in ("ID", "Parent") and key not in drop:
            out.append(kv)
    return ";".join(out)


def join_attrs(*parts):
    return ";".join(x for x in parts if x)


def span(block):
    cds = [c for c in block[2] if c[2] == "CDS"]
    return block[0][0], block[0][6], min(int(c[3]) for c in cds), max(int(c[4]) for c in cds)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--baskets", required=True)
    ap.add_argument("--overlap-votes", default="", help="basket_overlaps.py output (optional; old order)")
    ap.add_argument("--decisions", required=True)
    ap.add_argument("--winners", required=True, help="winners.tsv of baskets_v2.py (support and presence of each panel gene)")
    ap.add_argument("--own-invalid", required=True, help="own_invalid_orf.tsv of decide_models.py")
    ap.add_argument("--mode", choices=["keep_similar", "prefer_panel_model"], required=True)
    ap.add_argument("--native-dir", required=True)
    ap.add_argument("--transfer-dir", required=True)
    ap.add_argument("--domain-dir", required=True)
    ap.add_argument("--genomes", required=True, help="genome list, the reference first; its order gives the indexes")
    ap.add_argument("--old-reference-gff", default="", help="reference GFF before its curation (optional)")
    ap.add_argument("--curation-removed-list", default="",
                    help="IDs of the genes the curation deleted (one per line, IDs of --old-reference-gff); without it, "
                         "every gene of --old-reference-gff whose ID is not in the reference GFF counts as deleted")
    ap.add_argument("--array-models", default="", help="decide_models.py's array_models.tsv (tandem arrays)")
    ap.add_argument("--members", default="", help="baskets_v2.py's members.tsv: genes left without a model keep their own")
    ap.add_argument("--protected-baskets", default="",
                    help="baskets_v2.py's protected_baskets.tsv: their models win 'one place, one gene' and are tagged "
                         "protected=True")
    ap.add_argument("--id-prefix", required=True, help="prefix of the new gene IDs, e.g. Tg")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    genomes = [l.split()[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    ref = genomes[0]

    basket_of, winners = {}, []
    for r in csv.DictReader(open(a.baskets), delimiter="\t"):
        ws = [w for w in r["winners"].split(";") if w]   # a basket without a winner: ()
        winners += ws
        for m in r["members"].split(","):
            basket_of[m] = tuple(ws)
    wset = set(winners)
    prot_w = set()   # panel genes of the protected baskets (the shield)
    if a.protected_baskets and os.path.exists(a.protected_baskets):
        for r in csv.DictReader(open(a.protected_baskets), delimiter="\t"):
            prot_w.update(w for w in r["winners"].split(";") if w)
    lost_overlap = {r["gene"] for r in csv.DictReader(open(a.overlap_votes), delimiter="\t")
                    if int(r["would_drop_in"]) > 0} if a.overlap_votes else set()
    dom = {}
    for g in genomes:
        for r in csv.DictReader(open(os.path.join(a.domain_dir, f"{g}_gene_classification.tsv")), delimiter="\t"):
            dom[r["gene_id"]] = r["has_domain"] == "True"
    dec = defaultdict(dict)   # genome -> panel gene -> (decision, own genes, alternative)
    chosen_orf = {}
    for r in csv.DictReader(open(a.decisions), delimiter="\t"):
        dec[r["genome"]][r["panel_gene"]] = (r[f"decision_{a.mode}"], [x for x in r["own_genes"].split(",") if x],
                                             r["alternative_member"])
        # ORF status of the chosen model (valid / gap / end / invalid); not for keep_original (a valid own gene)
        if r.get("model_orf_prefer_panel_model") and r[f"decision_{a.mode}"] != "keep_original":
            chosen_orf[(r["genome"], r["panel_gene"])] = r["model_orf_prefer_panel_model"]

    info = {r["gene"]: (int(r["present_in"]), int(r["support"])) for r in csv.DictReader(open(a.winners), delimiter="\t")}
    own_invalid = {r["gene"] for r in csv.DictReader(open(a.own_invalid), delimiter="\t")}

    # native genes (all of them: blocks for the genes written, spans for the place check)
    native, owner, nspan = {}, {}, defaultdict(list)
    for g in genomes:
        path = next(iter((glob.glob(os.path.join(a.native_dir, g + ".gff3")) + glob.glob(os.path.join(a.native_dir, g + ".gff")))))
        allids = set()
        for line in open(path):
            f = line.split("\t")
            if len(f) > 8 and "Parent=" not in f[8]:
                m = ID_RE.search(f[8])
                if m:
                    allids.add(m.group(1))
        for k, b in gene_blocks(path, allids).items():
            native[k], owner[k] = b, g
            c, st, s, e = span(b)
            nspan[(g, c, st)].append((s, e, k))
    print(f"{len(native)} native genes", file=sys.stderr)

    # transfers needed: winners and alternative members, by source genome
    need = defaultdict(set)
    for g in genomes:
        for w, (d, own, alt) in dec[g].items():
            if d in ("replace", "fill", "target_invalid") and owner.get(w) != g:
                need[(owner.get(w), g)].add(w)
            elif d in ("alternative", "rescued") and alt and owner.get(alt) != g:
                need[(owner.get(alt), g)].add(alt)
    # tandem arrays (decide_models.py's array_models.tsv): each genome's own copies, and the copies added by transfer
    arr, arr_key = defaultdict(lambda: defaultdict(list)), {}
    if a.array_models and os.path.exists(a.array_models):
        for r in csv.DictReader(open(a.array_models), delimiter="\t"):
            arr[r["winners"]][r["genome"]].append((r["member"], r["model"], r["orf_status"]))
            if r["model"] == "transfer" and owner.get(r["member"]) != r["genome"]:
                need[(owner.get(r["member"]), r["genome"])].add(r["member"])
        for key in arr:
            for x in key.split(";"):
                arr_key[x] = key
    need[(None, None)] = set()
    tr = {}
    for (src, tgt), ids in need.items():
        if not src or not ids:
            continue
        path = os.path.join(a.transfer_dir, f"{src}_onto_{tgt}_liftoff.gff3_polished")
        for k, b in gene_blocks(path, ids, primary_only=True).items():
            tr[(k, tgt)] = b
    # curation: genes of the old reference GFF that the curation removed
    removed_cur = defaultdict(list)
    if a.old_reference_gff:
        cur_ids = {k for k in native if owner[k] == ref}
        # the curators' own list of deleted genes, when given: only these count (the old GFF gives their places);
        # without it, every old gene whose ID is not in the reference GFF counts (renamed genes would count too)
        listed = ({l.split()[0] for l in open(a.curation_removed_list) if l.strip() and not l.startswith("#")}
                  if a.curation_removed_list else None)
        old = {}
        for line in open(a.old_reference_gff):
            f = line.split("\t")
            if len(f) > 8 and f[2] in ("gene", "protein_coding_gene") and "Parent=" not in f[8]:
                m = ID_RE.search(f[8])
                if m and (m.group(1) in listed if listed is not None else m.group(1) not in cur_ids):
                    removed_cur[(f[0], f[6])].append((int(f[3]), int(f[4])))
        ref_models = {}
        want = defaultdict(set)
        for w in winners:
            if owner.get(w) and owner[w] != ref:
                want[owner[w]].add(w)
        for src, ids in want.items():
            path = os.path.join(a.transfer_dir, f"{src}_onto_{ref}_liftoff.gff3_polished")
            for k, b in gene_blocks(path, ids, primary_only=True).items():
                ref_models[k] = span(b)
    lost_curation = set()
    for w, (c, st, s, e) in (ref_models.items() if a.old_reference_gff else []):
        for s1, e1 in removed_cur.get((c, st), []):
            if min(e, e1) - max(s, s1) + 1 >= 0.5 * (e - s + 1):
                lost_curation.add(w)
                break

    def is_seleno(w):
        b = native.get(w)
        m = DESC_RE.search(b[0][8]) if b else None
        return bool(m and "selenoprotein" in m.group(1).lower())

    # selenoproteins that fall back to their own transfer need it loaded
    extra = defaultdict(set)
    for g in genomes:
        for w, (d, own, alt) in dec[g].items():
            if d in ("alternative", "rescued") and is_seleno(w) and owner.get(w) != g and (w, g) not in tr:
                extra[(owner[w], g)].add(w)
    for (src, tgt), ids in extra.items():
        path = os.path.join(a.transfer_dir, f"{src}_onto_{tgt}_liftoff.gff3_polished")
        for k, b in gene_blocks(path, ids, primary_only=True).items():
            tr[(k, tgt)] = b

    gone = lost_overlap | lost_curation   # panel genes removed from every genome

    # choose what each genome writes
    out = defaultdict(list)   # genome -> (contig, start, block, attrs)
    stats = defaultdict(Counter)
    removed_rows = []
    for w in winners:
        if w in lost_overlap:
            removed_rows.append([w, "overlap"])
            continue
        if w in lost_curation:
            removed_rows.append([w, "curation-removed reference gene"])
            continue
        for g in genomes:
            d, own, alt = dec[g].get(w, ("absent", [], ""))
            block, origin, orig_id, frm = None, None, "", w
            if d == "keep_original" and own:
                block, origin, orig_id, frm = native.get(own[0]), "original", own[0], own[0]
            elif d == "keep_divergent" and own:
                big = max(own, key=lambda k: span(native[k])[3] - span(native[k])[2] if k in native else -1)
                block, origin, orig_id, frm = native.get(big), "original_divergent", ",".join(own), big
            elif d == "target (own model)":
                block, origin, orig_id = native.get(w), "original", w
            elif d in ("replace", "fill", "target_invalid"):
                block = native.get(w) if owner.get(w) == g else tr.get((w, g))
                origin = {"replace": "replaced", "fill": "filled", "target_invalid": "target_invalid_orf"}[d]
                orig_id = ",".join(own) if d == "replace" else ""
            elif d in ("alternative", "rescued") and alt and is_seleno(w):
                # a selenoprotein reads its internal TGA as selenocysteine: keep the panel gene's own model
                block = tr.get((w, g))
                origin, orig_id = "target_invalid_orf", ",".join(own)
            elif d in ("alternative", "rescued") and alt:
                block = native.get(alt) if owner.get(alt) == g else tr.get((alt, g))
                origin, orig_id, frm = d, ",".join(own), alt
            if block is None:
                stats[g]["absent"] += 1
                continue
            if origin in ("replaced", "filled", "target_invalid_orf", "alternative", "rescued") and owner.get(frm) != g:
                c, st, s, e = span(block)
                clash = False
                for s1, e1, k in nspan.get((g, c, st), []):
                    other = basket_of.get(k)
                    if other is not None and all(x in gone for x in other):
                        continue   # a gene of a basket removed as false (or without a winner) does not hold its place
                    if other != basket_of.get(w) and min(e, e1) - max(s, s1) + 1 >= 0.5 * (e - s + 1):
                        clash = True
                        break
                if clash:
                    stats[g]["absent (transfer on another basket's gene)"] += 1
                    continue
            invalid = (origin == "target_invalid_orf" or (frm in own_invalid and owner.get(frm) == g)
                       or chosen_orf.get((g, w), "valid") in ("gap", "end", "invalid"))
            attrs = [f"original_ID={orig_id}" if orig_id else "", f"model_from={frm}", f"model_origin={origin}",
                     f"basket={w}", f"has_domain={dom.get(frm, False)}"]
            if owner.get(w) != ref:
                attrs.append("not_in_reference_annotation=True")
            if invalid:
                attrs.append("orf_invalid=True")
            cvm = re.search(r"(?:^|;)coverage=([\d.]+)", block[0][8])
            if origin in ("replaced", "filled", "target_invalid_orf", "alternative", "rescued") and cvm and float(cvm.group(1)) < 0.8:
                attrs.append(f"low_coverage_transfer=True;transfer_coverage={cvm.group(1)}")
            if is_seleno(w):
                attrs.append("selenoprotein=True")
            dm = DESC_RE.search(block[0][8])
            if dm:
                attrs.append(f"description={dm.group(1)}")
            pri = (owner.get(w) == ref, not invalid, info.get(w, (0, 0))[0], info.get(w, (0, 0))[1])
            out[g].append((block[0][0], int(block[0][3]), block, ";".join(x for x in attrs if x), pri, w, origin))
            stats[g][origin] += 1

    # tandem arrays: the genome's own confirmed copies and the valid transfers of other genomes' copies are added where
    # they do not overlap (same strand) a model already written for the basket, nor each other. Exception: a genome's
    # own copy with a valid ORF replaces the array's models over it when every one of them has a broken ORF (not a
    # selenoprotein); the replaced models are named in replaced_broken_model=.
    for key, per_g in arr.items():
        ws = set(key.split(";"))
        w0 = key.split(";")[0]
        if ws & gone:
            continue
        for g, extra in per_g.items():
            placed = [(span(row[2]), row) for row in out[g] if row[5] in ws]
            used = {re.search(r"model_from=([^;]+)", row[3]).group(1) for row in out[g] if row[5] in ws}
            for m, kind, status in extra:
                if m in used:
                    continue
                block = native.get(m) if kind == "own" else tr.get((m, g))
                if block is None:
                    continue
                sp = span(block)
                hits = [row for q, row in placed
                        if sp[0] == q[0] and sp[1] == q[1] and min(sp[3], q[3]) - max(sp[2], q[2]) + 1 > 0]
                replaced = []
                if hits:
                    if (kind == "own" and status == "valid"
                            and all("orf_invalid=True" in row[3] and "selenoprotein=True" not in row[3] for row in hits)):
                        replaced = hits
                    else:
                        continue
                if kind != "own":
                    c, st, s0, e0 = sp
                    if any(basket_of.get(k) != basket_of.get(w0) and min(e0, e1) - max(s0, s1) + 1 >= 0.5 * (e0 - s0 + 1)
                           for s1, e1, k in nspan.get((g, c, st), [])
                           if not (basket_of.get(k) is not None and all(x in gone for x in basket_of[k]))):
                        continue
                invalid = status in ("gap", "end", "invalid")
                origin = "original" if kind == "own" else "filled"
                attrs = [f"original_ID={m}" if kind == "own" else "", f"model_from={m}", f"model_origin={origin}",
                         f"basket={w0}", f"has_domain={dom.get(m, False)}", "tandem_array_copy=True"]
                if owner.get(w0) != ref:
                    attrs.append("not_in_reference_annotation=True")
                if invalid:
                    attrs.append("orf_invalid=True")
                if replaced:
                    attrs.append("replaced_broken_model=" + ",".join(
                        re.search(r"model_from=([^;]+)", row[3]).group(1) for row in replaced))
                    out[g] = [row for row in out[g] if not any(row is r for r in replaced)]
                    placed = [(q, row) for q, row in placed if not any(row is r for r in replaced)]
                    stats[g]["tandem array copy replaced a broken model"] += len(replaced)
                dm = DESC_RE.search(block[0][8])
                if dm:
                    attrs.append(f"description={dm.group(1)}")
                pri = (owner.get(w0) == ref, not invalid, info.get(w0, (0, 0))[0], info.get(w0, (0, 0))[1])
                row = (block[0][0], int(block[0][3]), block, ";".join(x for x in attrs if x), pri, w0, origin)
                out[g].append(row)
                placed.append((sp, row))
                used.add(m)
                stats[g]["tandem array copy (" + ("own" if kind == "own" else "transfer") + ")"] += 1

    # one place, one gene: two models on the same strand overlapping >= 50% of the smaller -> keep the one of higher
    # priority (reference basket, valid ORF, present in more genomes, more support)
    stacked_rows, dropped = [], defaultdict(list)
    for g in genomes:
        keep, by = [], defaultdict(list)
        # reference baskets first (the reference is never edited), then protected baskets, then the usual priority
        for row in sorted(out[g], key=lambda x: (x[4][0], x[5] in prot_w, x[4]), reverse=True):
            c, st, s0, e0 = span(row[2])
            clash = None
            for s1, e1, w1 in by[(c, st)]:
                if min(e0, e1) - max(s0, s1) + 1 >= 0.5 * min(e0 - s0 + 1, e1 - s1 + 1):
                    clash = w1
                    break
            if clash:
                stacked_rows.append([g, row[5], row[6], clash])
                dropped[g].append(row)
                stats[g][row[6]] -= 1
                stats[g]["dropped: same place as another panel gene"] += 1
                continue
            by[(c, st)].append((s0, e0, row[5]))
            keep.append(row)
        out[g] = keep

    # genes left without a model: each member of a basket with a winner that no model of its basket covers keeps its
    # own model when it passes the quality check and does not take the place of a model already written
    left_rows = []
    if a.members:
        check, rank = {}, {}
        for r in csv.DictReader(open(a.members), delimiter="\t"):
            check[r["gene"]] = r["keep_check"]
            rank[r["gene"]] = (float(r["psauron"]) if r["psauron"] else -1.0, int(r["protein_aa"] or 0))
        members_of = defaultdict(list)
        for m, ws in basket_of.items():
            if ws and owner.get(m):
                members_of[(owner[m], ws)].append(m)
        for (g, ws), ms in sorted(members_of.items()):
            if set(ws) & gone:
                continue   # the panel gene was removed (removed_panel_genes.tsv)
            # best PSAURON first, then the largest protein: in a basket without a model, the first one that passes
            # comes back and the others are pieces of it
            for m in sorted(ms, key=lambda m: (-rank.get(m, (-1.0, 0))[0], -rank.get(m, (-1.0, 0))[1], m)):
                if m not in native:
                    continue
                sp = span(native[m])
                rows_b = [row for row in out[g] if row[5] in ws]
                if any(sp[0] == q[0] and min(sp[3], q[3]) - max(sp[2], q[2]) + 1 > 0
                       for q in (span(row[2]) for row in rows_b)):
                    continue   # covered by a model of its basket
                w0 = next((w for w in ws if m in dec[g].get(w, ("", [], ""))[1]), ws[0])
                if rows_b:   # its basket has a model in this genome: a piece of that gene (copies left in baskets_v2.py)
                    q = span(rows_b[0][2])
                    left_rows.append([g, m, w0, f"piece of a gene that has its model here ({q[0]}:{q[2]}-{q[3]})"])
                    continue
                why = check.get(m, "not in members.tsv")
                if not why:
                    clash = next((row for row in out[g] for q in [span(row[2])]
                                  if q[0] == sp[0] and q[1] == sp[1] and
                                  min(sp[3], q[3]) - max(sp[2], q[2]) + 1 >= 0.5 * min(sp[3] - sp[2] + 1, q[3] - q[2] + 1)),
                                 None)
                    if clash:
                        why = f"same place as the model of {clash[5]}"
                if why:
                    left_rows.append([g, m, w0, why])
                    continue
                block = native[m]
                invalid = m in own_invalid
                attrs = [f"original_ID={m}", f"model_from={m}", "model_origin=original_other_place", f"basket={w0}",
                         f"has_domain={dom.get(m, False)}"]
                if owner.get(w0) != ref:
                    attrs.append("not_in_reference_annotation=True")
                if invalid:
                    attrs.append("orf_invalid=True")
                dm = DESC_RE.search(block[0][8])
                if dm:
                    attrs.append(f"description={dm.group(1)}")
                pri = (owner.get(w0) == ref, not invalid, info.get(w0, (0, 0))[0], info.get(w0, (0, 0))[1])
                out[g].append((block[0][0], int(block[0][3]), block, ";".join(attrs), pri, w0, "original_other_place"))
                stats[g]["original_other_place"] += 1
        for g, m, w0, why in left_rows:
            stats[g]["left out: " + ("PSAURON" if why.startswith("PSAURON") else "same place as another model"
                                     if why.startswith("same place") else "piece of a gene that has its model here"
                                     if why.startswith("piece") else why)] += 1

    def write_gff(path, rows, tag):
        with open(path, "w") as fh:
            fh.write("##gff-version 3\n")
            for n, (c, s, (gf, tf, parts), attrs, *_rest) in enumerate(sorted(rows, key=lambda x: (x[0], x[1])), 1):
                if len(_rest) > 1 and _rest[1] in prot_w:
                    attrs += ";protected=True"
                nid = f"{a.id_prefix}{gi}{tag}{n:05d}"
                tid = f"{nid}.t1"
                # the new tags first, then every attribute of the source line (only ID/Parent, and keys set here, replaced)
                set_keys = {kv.split("=", 1)[0] for kv in attrs.split(";") if "=" in kv}
                fh.write("\t".join(gf[:2] + ["gene"] + gf[3:8] + [join_attrs(f"ID={nid};{attrs}", keep_attrs(gf[8], set_keys))]) + "\n")
                ttype = "transcript" if "pseudogene=True" in gf[8] else "mRNA"   # a pseudogene: gene + transcript + exon
                fh.write("\t".join(tf[:2] + [ttype] + tf[3:8] + [join_attrs(f"ID={tid};Parent={nid}", keep_attrs(tf[8]))]) + "\n")
                k = Counter()
                for p in sorted(parts, key=lambda x: (x[2] != "exon", int(x[3]))):
                    k[p[2]] += 1
                    fh.write("\t".join(p[:8] + [join_attrs(f"ID={tid}.{p[2].lower()}{k[p[2]]};Parent={tid}", keep_attrs(p[8]))]) + "\n")

    # write
    for gi, g in enumerate(genomes, 1):
        write_gff(os.path.join(a.out_dir, f"{g}_basket_annotation.gff"), out[g], "g")
        write_gff(os.path.join(a.out_dir, f"{g}_stacked_dropped.gff"), dropped[g], "s")
    with open(os.path.join(a.out_dir, "stacked_dropped.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "panel_gene_dropped", "model_origin", "kept_panel_gene"])
        w.writerows(stacked_rows)
    with open(os.path.join(a.out_dir, "left_out.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "gene", "panel_gene", "reason"])
        w.writerows(left_rows)
    with open(os.path.join(a.out_dir, "removed_panel_genes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["panel_gene", "reason"])
        w.writerows(removed_rows)
    keys = sorted({k for c in stats.values() for k in c})
    with open(os.path.join(a.out_dir, "summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["genome", "genes"] + keys)
        for g in genomes:
            w.writerow([g, sum(v for k, v in stats[g].items() if not k.startswith(("absent", "dropped", "left out")))] + [stats[g][k] for k in keys])
    print(f"mode {a.mode}: panel genes {len(winners)}, removed for overlap {len(lost_overlap & wset)}, "
          f"removed as curation-removed {len(lost_curation - lost_overlap)}")
    for g in genomes:
        n = sum(v for k, v in stats[g].items() if not k.startswith(("absent", "dropped", "left out")))
        print(f"  {g:10} {n:6} genes  " + "  ".join(f"{k}: {v}" for k, v in sorted(stats[g].items())))


if __name__ == "__main__":
    main()
