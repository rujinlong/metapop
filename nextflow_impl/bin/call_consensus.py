#!/usr/bin/env python3
"""
Call consensus bases from multiple SNP count files.
Adapted from metapop_snp_call.py choose_consensus() function.
"""

import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(description="Call consensus bases from SNP counts.")
    parser.add_argument("--input_dir", required=True, help="Directory containing SNP count files")
    parser.add_argument("--output", required=True, help="Output consensus file")
    
    args = parser.parse_args()
    
    # Find all SNP count files
    snp_files = glob.glob(os.path.join(args.input_dir, "*_snp_counts.tsv"))
    
    # Aggregate counts across all samples
    consensus = {}  # (contig, pos) -> [A, T, C, G]
    
    for snp_file in snp_files:
        with open(snp_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                contig = parts[0]
                pos = int(parts[1])
                counts = [int(x) for x in parts[2:6]]
                
                key = (contig, pos)
                if key not in consensus:
                    consensus[key] = [0, 0, 0, 0]
                for i in range(4):
                    consensus[key][i] += counts[i]
    
    # Determine consensus base
    bases = ['A', 'T', 'C', 'G']
    
    with open(args.output, 'w') as out:
        print("contig\tpos\tref_base\tA\tT\tC\tG", file=out)
        
        for (contig, pos), counts in sorted(consensus.items()):
            max_idx = counts.index(max(counts))
            ref_base = bases[max_idx]
            print(f"{contig}\t{pos}\t{ref_base}\t{counts[0]}\t{counts[1]}\t{counts[2]}\t{counts[3]}", file=out)

if __name__ == "__main__":
    main()
