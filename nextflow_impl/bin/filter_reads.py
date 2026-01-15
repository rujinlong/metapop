#!/usr/bin/env python3

import pysam
import argparse
import re
import sys
import os

def parse_entry(entry, out_handle, is_global, min_length, min_pct_id):
    if not entry.has_tag("MD"):
        # If no MD tag, we can't calculate %ID easily without ref. 
        # Metapop logic implies MD tags are present (calmd run before).
        return None
    
    # Get the match/mismatch info for the read
    mdz_seg = entry.get_tag("MD")
    match_count = re.findall('[0-9]+', mdz_seg)
    
    if is_global:
        # pct ID on global is %match/read length
        # Note: Metapop original used total_count = entry.query_length
        # But wait, original code: total_count = entry.query_length
        sum_matches = 0
        for num in match_count:
            sum_matches += int(num)
        
        total_count = entry.query_length
        
        if total_count == 0: # Avoid div by zero
            return None
            
        pct_id = (sum_matches / total_count) * 100
        
        if pct_id >= min_pct_id and total_count >= min_length:
            out_handle.write(entry)
            
    else:
        # pct ID on local is %match/alignment length
        sum_matches = 0
        for num in match_count:
            sum_matches += int(num)
            
        # Total count here is num matches + num mismatches?
        # Metapop original: len(''.join([i for i in mdz_seg if not i.isdigit()])) + sum
        # The non-digits in MD tag are deleted bases from ref (A, T, C, G, ^).
        # Actually this calculates matches + deletions + mismatches? 
        # Alignment length usually includes insertions too (CIGAR).
        # Let's stick strictly to Metapop's logic to maintain identical behavior.
        
        non_digit_len = len(''.join([i for i in mdz_seg if not i.isdigit()]))
        # Note: MD tag doesn't show insertions (I). 
        # If Metapop uses this formula, we use it too.
        total_count = non_digit_len + sum_matches
        
        if total_count == 0:
            return None
            
        pct_id = (sum_matches / total_count) * 100
        
        # Original: if pct_id >= min_pct_id and total_count >= min_length:
        if pct_id >= min_pct_id and total_count >= min_length:
            out_handle.write(entry)

def main():
    parser = argparse.ArgumentParser(description="Filter BAM by %ID and length.")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument("--min_len", type=int, default=50, help="Minimum read length")
    parser.add_argument("--min_pct_id", type=float, default=95, help="Minimum percent identity")
    parser.add_argument("--global_id", action="store_true", help="Use global percent identity (matches/read_len)")
    
    args = parser.parse_args()
    
    # Open BAMs
    # 'rb' for read bam, 'wb' for write bam
    input_reads = pysam.AlignmentFile(args.input, "rb")
    output_reads = pysam.AlignmentFile(args.output, "wb", template=input_reads)
    
    count_pass = 0
    count_total = 0
    
    for read in input_reads:
        count_total += 1
        # Check if unmapped
        if read.is_unmapped:
            continue
            
        parse_entry(read, output_reads, args.global_id, args.min_len, args.min_pct_id)
        
    input_reads.close()
    output_reads.close()
    
    # We might want to print stats to stderr
    # print(f"Processed {count_total} reads.", file=sys.stderr)

if __name__ == "__main__":
    main()
