#!/usr/bin/env python3

import sys
import math
import argparse
import csv

def calculate_stats(counts_dict, length, percentage_truncation):
    # logic from breadth_and_depth_output in metapop_filter.py
    
    num_pos_covered = 0
    for depth in counts_dict:
        num_pos_covered += counts_dict[depth]
        
    breadth = (num_pos_covered / length) * 100
    
    # Calculate TAD (Truncated Average Depth)
    cutoff = math.floor(length / percentage_truncation)
    
    if num_pos_covered <= cutoff:
        TAD = 0.0
    else:
        # Add zeroes for uncovered positions
        count_zeroes = length - num_pos_covered
        # counts_dict is passed by reference, but we are about to modify it destructively for TAD calc.
        # Let's make a local copy or be careful.
        # Actually in the original code, it modifies the dict. Since we process one contig at a time and reset, it's fine.
        counts_dict[0] = count_zeroes
        
        observed_depths = list(counts_dict.keys())
        observed_depths.sort()
        
        # Remove lower values
        depths_removed_so_far = 0
        working_counts = counts_dict.copy() # Work on copy to be safe
        
        for depth in observed_depths:
            obs = working_counts[depth]
            if depths_removed_so_far + obs >= cutoff:
                working_counts[depth] -= (cutoff - depths_removed_so_far)
                break
            else:
                depths_removed_so_far += obs
                working_counts[depth] = 0
                
        # Remove upper values
        depths_removed_so_far = 0
        observed_depths.sort(reverse=True)
        for depth in observed_depths:
            obs = working_counts[depth]
            if depths_removed_so_far + obs >= cutoff:
                working_counts[depth] -= (cutoff - depths_removed_so_far)
                break
            else:
                depths_removed_so_far += obs
                working_counts[depth] = 0
                
        # Calc TAD
        tad_sum = 0
        for depth in working_counts:
            tad_sum += (depth * working_counts[depth])
            
        total_count = length - (2 * cutoff)
        
        if total_count > 0:
            TAD = round(tad_sum / total_count, 6)
        else:
            TAD = 0.0
            
    return breadth, TAD

def main():
    parser = argparse.ArgumentParser(description="Calculate coverage stats from samtools depth.")
    parser.add_argument("--lengths", required=True, help="File containing contig lengths (two columns: contig length)")
    parser.add_argument("--min_cov", type=float, default=20.0, help="Min breadth coverage %%")
    parser.add_argument("--min_dep", type=float, default=10.0, help="Min TAD")
    parser.add_argument("--trunc", type=float, default=10.0, help="Truncation percentile (e.g. 10)")
    parser.add_argument("--output_stats", required=True, help="Output file for stats")
    parser.add_argument("--output_pass", required=True, help="Output file for passing genomes list")
    
    args = parser.parse_args()
    
    # Load lengths
    contig_lengths = {}
    with open(args.lengths, 'r') as f:
        # Fasta index (.fai) or just tab delimited
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                contig_lengths[parts[0]] = int(parts[1])
    
    # Process stdin (samtools depth output)
    # Format: contig pos depth
    
    current_contig = None
    current_depths_counts = {} # depth -> count
    
    stats_fh = open(args.output_stats, 'w')
    pass_fh = open(args.output_pass, 'w')
    
    # Write header for stats
    print("contig\tnum_pos_covered\tbreadth\tTAD", file=stats_fh)
    
    for line in sys.stdin:
        parts = line.strip().split('\t')
        if len(parts) < 3: 
            continue
            
        contig = parts[0]
        # pos = int(parts[1]) # not needed for aggregation
        depth = int(parts[2])
        
        if contig != current_contig:
            if current_contig is not None:
                # Process previous contig
                if current_contig in contig_lengths:
                    length = contig_lengths[current_contig]
                    breadth, tad = calculate_stats(current_depths_counts, length, args.trunc)
                    print(f"{current_contig}\t{sum(current_depths_counts.values())}\t{breadth}\t{tad}", file=stats_fh)
                    
                    if breadth >= args.min_cov and tad >= args.min_dep:
                        print(current_contig, file=pass_fh)
                else:
                    # Contig in depth but not in lengths file?
                    sys.stderr.write(f"Warning: Contig {current_contig} found in depth but not in lengths file.\n")
                    
            current_contig = contig
            current_depths_counts = {}
            
        current_depths_counts[depth] = current_depths_counts.get(depth, 0) + 1
        
    # Process last contig
    if current_contig is not None and current_contig in contig_lengths:
        length = contig_lengths[current_contig]
        breadth, tad = calculate_stats(current_depths_counts, length, args.trunc)
        print(f"{current_contig}\t{sum(current_depths_counts.values())}\t{breadth}\t{tad}", file=stats_fh)
        
        if breadth >= args.min_cov and tad >= args.min_dep:
            print(current_contig, file=pass_fh)
            
    stats_fh.close()
    pass_fh.close()

if __name__ == "__main__":
    main()
