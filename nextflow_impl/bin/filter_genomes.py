#!/usr/bin/env python3

import pysam
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Filter BAM to retain only specific genomes.")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument("--genomes", required=True, help="File containing list of allowed genomes (one per line)")
    
    args = parser.parse_args()
    
    # Load allowed genomes
    allowed = set()
    with open(args.genomes, 'r') as f:
        for line in f:
            allowed.add(line.strip().split()[0])
            
    input_reads = pysam.AlignmentFile(args.input, "rb")
    output_reads = pysam.AlignmentFile(args.output, "wb", template=input_reads)
    
    count = 0
    kept = 0
    
    for read in input_reads:
        count += 1
        # get_reference_name returns string or None
        ref = read.reference_name
        
        if ref and ref in allowed:
            output_reads.write(read)
            kept += 1
            
    input_reads.close()
    output_reads.close()
    
    # print(f"Processed {count} reads, kept {kept}.", file=sys.stderr)

if __name__ == "__main__":
    main()
