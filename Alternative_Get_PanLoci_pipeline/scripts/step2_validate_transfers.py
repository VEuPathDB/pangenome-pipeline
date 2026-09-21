#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import csv
import glob
import tempfile
from collections import defaultdict

def load_fasta(fasta_path):
    """
    FASTA parser that loads the entire genome into memory.
    Returns a dictionary of {header: sequence_string}.
    """
    fasta_dict = {}
    if not os.path.exists(fasta_path):
        print(f"Error: FASTA file not found at '{fasta_path}'", file=sys.stderr)
        return None
        
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    fasta_dict[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]  # First word of header
                current_seq = []
            else:
                current_seq.append(line.upper())
                
        if current_header:
            fasta_dict[current_header] = "".join(current_seq)
            
    return fasta_dict

def extract_sequence(fasta_dict, chrom, start, end, strand="+"):
    """Extracts a sequence substring from the in-memory genome FASTA dictionary."""
    if not fasta_dict or chrom not in fasta_dict:
        return None
        
    # Ensure coordinates are within bounds
    chrom_len = len(fasta_dict[chrom])
    start = max(0, min(start, chrom_len))
    end = max(0, min(end, chrom_len))
    
    if start >= end:
        return None
        
    seq = fasta_dict[chrom][start:end]
    
    if strand == "-":
        # Reverse complement (efficient translation)
        comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        seq = seq.translate(comp)[::-1]
        
    return seq

def run_batch_minimap2(queries, targets, output_prefix, minimap_bin="minimap2", cactus_img=None):
    """
    Executes a database-wide batch minimap2 alignment using permissive parameters
    optimized for AT-rich genomes (-c -N 100 -p 0.1).
    """
    q_file = f"{output_prefix}_queries.fa"
    t_file = f"{output_prefix}_targets.fa"
    
    # Write multi-FASTA for queries
    with open(q_file, 'w') as qf:
        for gene_id, seq in queries.items():
            if seq:
                qf.write(f">{gene_id}\n{seq}\n")
                
    # Write multi-FASTA for targets
    with open(t_file, 'w') as tf:
        for gene_id, seq in targets.items():
            if seq:
                tf.write(f">{gene_id}\n{seq}\n")
                
    # Permissive batch parameters reduce early heuristic filtering in AT-rich genomes.
    if cactus_img:
        cmd = [
            "apptainer", "exec", cactus_img, "minimap2",
            "-c", "-N", "100", "-p", "0.1", t_file, q_file
        ]
    else:
        cmd = [minimap_bin, "-c", "-N", "100", "-p", "0.1", t_file, q_file]
        
    paf_lines = []
    try:
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        paf_output = result.stdout.strip()
        if paf_output:
            paf_lines = paf_output.split('\n')
    except subprocess.CalledProcessError as e:
        print(f"[Warning] Batch minimap2 failed: {e.stderr}", file=sys.stderr)
    finally:
        # Clean up batch FASTA files immediately
        if os.path.exists(q_file):
            os.remove(q_file)
        if os.path.exists(t_file):
            os.remove(t_file)
            
    # Parse batch results (keep the alignment with highest matches for each orthologous pair)
    results = {}
    for line in paf_lines:
        if not line.strip():
            continue
        cols = line.split('\t')
        if len(cols) < 11:
            continue
            
        q_name = cols[0]
        t_name = cols[5]
        
        # Keep only orthologous self-locus alignments (where query ID matches target ID)
        if q_name == t_name:
            query_len = int(cols[1])
            q_start = int(cols[2])
            q_end = int(cols[3])
            matches = int(cols[9])
            block_len = int(cols[10])
            
            coverage = round((q_end - q_start) / query_len, 4) if query_len > 0 else 0.0
            identity = round((matches / block_len) * 100, 2) if block_len > 0 else 0.0
            
            if q_name not in results or matches > results[q_name]["matches"]:
                results[q_name] = {
                    "query_len": query_len,
                    "target_len": int(cols[6]),
                    "query_coverage": coverage,
                    "percent_identity": identity,
                    "alignment_length": block_len,
                    "matches": matches
                }
            
    return results

def run_minimap2_pairwise_fallback(seq_ref, seq_target, gene_id, minimap_bin="minimap2", cactus_img=None):
    """
    Executes an isolated 1-to-1 pairwise fallback alignment for locus pairs
    not recovered by the batch alignment.
    """
    # Create thread-safe uniquely named temporary FASTA files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f_ref, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f_target:
        f_ref.write(f">{gene_id}_ref\n{seq_ref}\n")
        f_target.write(f">{gene_id}_target\n{seq_target}\n")
        f_ref_name = f_ref.name
        f_target_name = f_target.name
        
    paf_lines = []
    try:
        # Use standard minimap2 parameters for isolated fallback alignments.
        if cactus_img:
            cmd = ["apptainer", "exec", cactus_img, "minimap2", "-c", f_target_name, f_ref_name]
        else:
            cmd = [minimap_bin, "-c", f_target_name, f_ref_name]
            
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        paf_output = res.stdout.strip()
        if paf_output:
            paf_lines = paf_output.split('\n')
                
    except Exception as e:
        print(f"[Warning] Pairwise fallback failed for {gene_id}: {e}", file=sys.stderr)
    finally:
        # Guarantee cleanup of temporary files
        if os.path.exists(f_ref_name):
            os.remove(f_ref_name)
        if os.path.exists(f_target_name):
            os.remove(f_target_name)
            
    for line in paf_lines:
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) < 11:
            continue
            
        query_len = int(fields[1])
        q_start = int(fields[2])
        q_end = int(fields[3])
        matches = int(fields[9])
        block_len = int(fields[10])
        
        coverage = round((q_end - q_start) / query_len, 4) if query_len > 0 else 0.0
        identity = round((matches / block_len) * 100, 2) if block_len > 0 else 0.0
        
        return {
            "query_len": query_len,
            "target_len": int(fields[6]),
            "query_coverage": coverage,
            "percent_identity": identity,
            "alignment_length": block_len,
            "matches": matches
        }
    return None

def classify_occupancy(
    coverage,
    identity,
    present_cov_cutoff=0.50,
    present_id_cutoff=50.0,
    partial_cov_cutoff=0.30,
    partial_id_cutoff=30.0,
):
    """
    Classifies locus occupancy based on sequence validation metrics:
      - PRESENT: coverage >= present_cov_cutoff and identity >= present_id_cutoff
      - PARTIAL: coverage >= partial_cov_cutoff and identity >= partial_id_cutoff
      - ABSENT: coverage < partial_cov_cutoff or identity < partial_id_cutoff
    """
    if coverage >= present_cov_cutoff and identity >= present_id_cutoff:
        return "PRESENT"
    elif coverage >= partial_cov_cutoff and identity >= partial_id_cutoff:
        return "PARTIAL"
    else:
        return "ABSENT"

def parse_samples(samples_tsv):
    """Parses samples.tsv to extract strain information."""
    samples = []
    if not os.path.exists(samples_tsv):
        print(f"Error: Metadata file not found at '{samples_tsv}'", file=sys.stderr)
        sys.exit(1)
        
    with open(samples_tsv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            samples.append({
                "strain_id": row["strain_id"].strip(),
                "genome_fasta": row["genome_fasta"].strip(),
                "annotation_gff": row["annotation_gff"].strip(),
                "analysis_order": row["analysis_order"].strip()
            })
    return samples

def main():
    parser = argparse.ArgumentParser(description="Step 2: validate transferred panloci with minimap2")
    parser.add_argument("-s", "--samples", default="samples.tsv", help="Path to samples.tsv metadata (default: samples.tsv)")
    parser.add_argument("-b", "--beddir", default="./candidate_transfers", help="Directory containing merged BED files (default: ./candidate_transfers)")
    parser.add_argument("-o", "--outdir", default="./validation", help="Output directory for validation reports (default: ./validation)")
    parser.add_argument("-m", "--minimap-bin", default="minimap2", help="Path to minimap2 executable (default: minimap2)")
    parser.add_argument("-c", "--cactus-img", help="Path to Cactus Apptainer SIF image (optional, runs minimap2 via Singularity)")
    parser.add_argument("--present-cov-cutoff", type=float, default=0.50, help="Coverage cutoff for PRESENT (default: 0.50)")
    parser.add_argument("--present-id-cutoff", type=float, default=50.0, help="Identity cutoff for PRESENT (default: 50.0)")
    parser.add_argument("--partial-cov-cutoff", type=float, default=0.30, help="Coverage cutoff for PARTIAL (default: 0.30)")
    parser.add_argument("--partial-id-cutoff", type=float, default=30.0, help="Identity cutoff for PARTIAL (default: 30.0)")
    
    args = parser.parse_args()
    
    print("======================================================================")
    print(" Step 2: Validate Transferred Panloci")
    print("======================================================================")
    print(f"  Metadata TSV:       {args.samples}")
    print(f"  BED Directory:      {args.beddir}")
    print(f"  Validation Output:  {args.outdir}")
    if args.cactus_img:
        print(f"  Cactus SIF Image:   {args.cactus_img}")
    else:
        print(f"  Minimap2 Bin:       {args.minimap_bin}")
    print(f"  PRESENT cutoffs:    Coverage >= {args.present_cov_cutoff}, Identity >= {args.present_id_cutoff}")
    print(f"  PARTIAL cutoffs:    Coverage >= {args.partial_cov_cutoff}, Identity >= {args.partial_id_cutoff}")
    print("======================================================================\n")
    
    # 1. Parse Samples
    samples = parse_samples(args.samples)
    print(f"[+] Loaded {len(samples)} strain records from metadata.")
    
    # 2. Locate and load reference seed strain (order = 1)
    ref_genome = next((s for s in samples if s["analysis_order"] == "1"), None)
    if not ref_genome:
        print("Error: Could not find reference genome with analysis_order = 1 in samples.tsv", file=sys.stderr)
        sys.exit(1)
        
    ref_id = ref_genome["strain_id"]
    ref_fasta_path = ref_genome["genome_fasta"]
    print(f"[+] Reference Seed Genome identified: {ref_id} ({ref_fasta_path})")
    
    print("[+] Loading reference genome FASTA into memory (cached once)...")
    ref_fasta = load_fasta(ref_fasta_path)
    if not ref_fasta:
        print(f"Error: Could not load reference FASTA from '{ref_fasta_path}'", file=sys.stderr)
        sys.exit(1)
    print(f"  --> Cached {len(ref_fasta)} reference chromosomes.")
    
    # Locate reference seed BED coordinates
    seed_beds = glob.glob(os.path.join(args.beddir, "*_seed.bed"))
    if not seed_beds:
        print(f"Error: Could not locate seed BED file (*_seed.bed) inside: {args.beddir}", file=sys.stderr)
        sys.exit(1)
    seed_bed_path = seed_beds[0]
    print(f"[+] Located seed reference BED file: {seed_bed_path}")
    
    # Parse seed coordinates
    ref_coords = {}
    with open(seed_bed_path, 'r', encoding='utf-8') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) >= 6:
                ref_coords[cols[3]] = {
                    "chrom": cols[0],
                    "start": int(cols[1]),
                    "end": int(cols[2]),
                    "strand": cols[5]
                }
    print(f"  --> Loaded {len(ref_coords)} reference coordinate definitions.")
    
    # 3. Iterate over target genomes and run pairwise validation
    os.makedirs(args.outdir, exist_ok=True)
    
    for s in samples:
        target_strain = s["strain_id"]
        if target_strain == ref_id:
            continue
            
        target_fasta_path = s["genome_fasta"]
        merged_bed_path = os.path.join(args.beddir, f"{ref_id}_to_{target_strain}_merged.bed")
        output_tsv_path = os.path.join(args.outdir, f"{ref_id}_to_{target_strain}.validation.tsv")
        
        print(f"\n[+] Validating coordinate transfers for: {ref_id} -> {target_strain}")
        
        if not os.path.exists(merged_bed_path):
            print(f"  [Warning] Merged BED file not found: {merged_bed_path}. Skipping target strain.", file=sys.stderr)
            continue
            
        # Parse target transfers
        target_transfers = []
        with open(merged_bed_path, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                if len(cols) >= 8:
                    target_transfers.append({
                        "chrom": cols[0],
                        "start": int(cols[1]),
                        "end": int(cols[2]),
                        "name": cols[3],
                        "strand": cols[5],
                        "orig_size": int(cols[6]),
                        "lifted_size": int(cols[7])
                    })
        
        print(f"  --> Loaded {len(target_transfers)} transferred locus intervals.")
        
        # Load target genome FASTA
        print(f"  --> Loading target genome FASTA: {target_fasta_path} ...")
        target_fasta = load_fasta(target_fasta_path)
        if not target_fasta:
            print(f"  [Warning] Failed to load FASTA for {target_strain}. Skipping.", file=sys.stderr)
            continue
            
        # 3.1 Extract sequences for all genes in batch
        queries = {}
        targets = {}
        for t in target_transfers:
            gene_id = t["name"]
            if gene_id not in ref_coords:
                continue
            
            r_coord = ref_coords[gene_id]
            queries[gene_id] = extract_sequence(ref_fasta, r_coord["chrom"], r_coord["start"], r_coord["end"], r_coord["strand"])
            targets[gene_id] = extract_sequence(target_fasta, t["chrom"], t["start"], t["end"], t["strand"])
            
        # Run one batch minimap2 alignment per target strain.
        print("  --> Running database-wide batch minimap2 alignment...")
        prefix_temp = os.path.join(args.outdir, f"temp_{ref_id}_to_{target_strain}")
        batch_results = run_batch_minimap2(queries, targets, prefix_temp, args.minimap_bin, args.cactus_img)
        
        # Run isolated fallback alignments for loci not recovered in the batch results.
        print("  --> Running isolated 1-to-1 pairwise fallback for missing loci...")
        fallback_results = {}
        fallback_count = 0
        
        for t in target_transfers:
            gene_id = t["name"]
            # If the locus is missing from batch results, try an isolated pairwise alignment.
            if gene_id not in batch_results:
                seq_ref = queries.get(gene_id)
                seq_target = targets.get(gene_id)
                if seq_ref and seq_target:
                    res_fb = run_minimap2_pairwise_fallback(seq_ref, seq_target, gene_id, args.minimap_bin, args.cactus_img)
                    if res_fb:
                        fallback_results[gene_id] = res_fb
                        fallback_count += 1
                        
        print(f"      [INFO] Recovered {fallback_count} loci via 1-to-1 pairwise fallback.")
        
        # Classify occupancy and prepare validation rows.
        target_transfers_dict = {t["name"]: t for t in target_transfers}
        validation_records = []
        success_count = 0
        partial_count = 0
        absent_count = 0
        
        for gene_id in ref_coords:
            ref_len = len(queries[gene_id]) if queries.get(gene_id) else 0
            
            t = target_transfers_dict.get(gene_id)
            if t:
                orig_size = t["orig_size"]
                lifted_size = t["lifted_size"]
                
                # Combine batch and fallback results
                res = None
                if gene_id in batch_results:
                    res = batch_results[gene_id]
                elif gene_id in fallback_results:
                    res = fallback_results[gene_id]
                    
                if res:
                    # Query Coverage is defined as the fraction of original donor CDS successfully transferred (lifted_size / orig_size)
                    cov = round(lifted_size / orig_size, 4) if orig_size > 0 else 0.0
                    ident = res["percent_identity"]
                    aln_len = res["alignment_length"]
                    t_len = res["target_len"]
                    status = classify_occupancy(
                        cov,
                        ident,
                        args.present_cov_cutoff,
                        args.present_id_cutoff,
                        args.partial_cov_cutoff,
                        args.partial_id_cutoff,
                    )
                else:
                    # No alignment was reported under any parameters
                    status = "ABSENT"
                    cov = 0.0
                    ident = 0.0
                    aln_len = 0
                    t_len = len(targets[gene_id]) if targets.get(gene_id) else 0
            else:
                # Completely missing from coordinate transfer -> ABSENT
                status = "ABSENT"
                cov = 0.0
                ident = 0.0
                aln_len = 0
                t_len = 0
                
            if status == "PRESENT":
                success_count += 1
            elif status == "PARTIAL":
                partial_count += 1
            else:
                absent_count += 1
                
            validation_records.append([
                gene_id,
                ref_len,
                t_len,
                aln_len,
                cov,
                ident,
                status
            ])
            
        # Write target validation output TSV
        with open(output_tsv_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["gene_id", "query_length", "target_length", "alignment_length", "query_coverage", "percent_identity", "validation_status"])
            writer.writerows(validation_records)
            
        print(f"      [SUCCESS] Written validation report: {output_tsv_path}")
        print(f"      Summary: PRESENT: {success_count} | PARTIAL: {partial_count} | ABSENT: {absent_count}")
        
        # Free target FASTA from memory to prevent leaks
        del target_fasta
        
    print("\n======================================================================")
    print(" Pangenome DNA sequence validation complete")
    print(f" Validation reports saved in: {args.outdir}")
    print("======================================================================\n")

if __name__ == "__main__":
    main()
