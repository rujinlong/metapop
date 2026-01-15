#!/usr/bin/env python3
"""
Refine SNP calls by counting alleles at each variant position.
Adapted from metapop_snp_call.py pileup() function.
"""

import pysam
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Refine SNP calls by counting alleles.")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--variants", required=True, help="Combined variants file (all samples)")
    parser.add_argument("--min_qual", type=int, default=20, help="Minimum base quality")
    parser.add_argument("--output", required=True, help="Output TSV file")
    
    args = parser.parse_args()
    
    # Load all variant positions
    variant_positions = {}  # contig -> set of positions
    with open(args.variants, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            contig = parts[0].split()[0]
            pos = int(parts[1])
            if contig not in variant_positions:
                variant_positions[contig] = set()
            variant_positions[contig].add(pos)
    
    # Open BAM file
    bam = pysam.AlignmentFile(args.bam, "rb")
    
    with open(args.output, 'w') as out:
        # Header
        print("contig\tpos\tA\tT\tC\tG", file=out)
        
        # Iterate through pileup at variant positions
        for contig in variant_positions:
            positions = sorted(variant_positions[contig])
            
            for pos in positions:
                # pysam uses 0-based coordinates
                counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
                
                try:
                    for pileupcolumn in bam.pileup(contig, pos-1, pos, 
                                                   min_base_quality=args.min_qual,
                                                   truncate=True):
                        if pileupcolumn.pos == pos - 1:  # 0-based
                            for pileupread in pileupcolumn.pileups:
                                if not pileupread.is_del and not pileupread.is_refskip:
                                    base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                                    if base in counts:
                                        counts[base] += 1
                except ValueError:
                    # Contig not in BAM
                    continue
                
                print(f"{contig}\t{pos}\t{counts['A']}\t{counts['T']}\t{counts['C']}\t{counts['G']}", file=out)
    
    bam.close()

if __name__ == "__main__":
    main()
