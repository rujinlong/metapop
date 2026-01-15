#!/usr/bin/env python3
"""
Associate SNPs with genes and calculate codon positions.
Simplified version adapted from metapop_snp_call.py associate_genes().
"""

import argparse
import os
import glob
import re
from collections import defaultdict

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def parse_prodigal_header(header):
    """Parse Prodigal FASTA header to extract gene info."""
    # Format: >contig_genenum # start # end # strand # metadata
    parts = header.split(' # ')
    if len(parts) < 4:
        return None
    
    gene_id = parts[0].lstrip('>')
    start = int(parts[1])
    end = int(parts[2])
    strand = int(parts[3])  # 1 or -1
    
    # Extract contig name (everything before last _)
    gene_parts = gene_id.rsplit('_', 1)
    contig = gene_parts[0] if len(gene_parts) > 1 else gene_id
    
    return {
        'gene_id': gene_id,
        'contig': contig,
        'start': start,
        'end': end,
        'strand': strand
    }

def load_genes(genes_file):
    """Load gene annotations from Prodigal output."""
    genes = defaultdict(list)  # contig -> list of gene dicts
    
    with open(genes_file, 'r') as f:
        current_gene = None
        current_seq = []
        
        for line in f:
            if line.startswith('>'):
                if current_gene and current_seq:
                    current_gene['sequence'] = ''.join(current_seq).upper()
                    genes[current_gene['contig']].append(current_gene)
                
                current_gene = parse_prodigal_header(line.strip())
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        # Last gene
        if current_gene and current_seq:
            current_gene['sequence'] = ''.join(current_seq).upper()
            genes[current_gene['contig']].append(current_gene)
    
    # Sort genes by start position
    for contig in genes:
        genes[contig].sort(key=lambda x: x['start'])
    
    return genes

def main():
    parser = argparse.ArgumentParser(description="Associate SNPs with genes.")
    parser.add_argument("--snp_dir", required=True, help="Directory containing SNP count files")
    parser.add_argument("--consensus", required=True, help="Consensus bases file")
    parser.add_argument("--genes", required=True, help="Prodigal genes file")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--min_obs", type=int, default=2, help="Min observations")
    parser.add_argument("--min_pct", type=float, default=1, help="Min percent")
    parser.add_argument("--output_genic", required=True, help="Genic SNPs output")
    parser.add_argument("--output_non_genic", required=True, help="Non-genic SNPs output")
    parser.add_argument("--output_genes", required=True, help="Corrected genes output")
    parser.add_argument("--output_codons", required=True, help="Codon counts output")
    
    args = parser.parse_args()
    
    # Load consensus bases
    consensus = {}  # (contig, pos) -> {'base': X, 'counts': [A,T,C,G]}
    with open(args.consensus, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            contig = parts[0]
            pos = int(parts[1])
            ref_base = parts[2]
            counts = [int(x) for x in parts[3:7]]
            consensus[(contig, pos)] = {'base': ref_base, 'counts': counts}
    
    # Load genes
    genes = load_genes(args.genes)
    
    # Process SNPs and associate with genes
    bases = ['A', 'T', 'C', 'G']
    min_pct_decimal = args.min_pct / 100.0
    
    genic_out = open(args.output_genic, 'w')
    non_genic_out = open(args.output_non_genic, 'w')
    
    # Headers
    print("contig\tpos\tref_base\talt_bases\tdepth\tgene_id\tstart\tend\tstrand\tcodon\tpos_in_codon\tsample", file=genic_out)
    print("contig\tpos\tref_base\talt_bases\tdepth\tsample", file=non_genic_out)
    
    # Process each SNP count file
    snp_files = glob.glob(os.path.join(args.snp_dir, "*_snp_counts.tsv"))
    
    for snp_file in snp_files:
        sample = os.path.basename(snp_file).replace("_snp_counts.tsv", "")
        
        with open(snp_file, 'r') as f:
            next(f)  # Skip header
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                
                contig = parts[0]
                pos = int(parts[1])
                counts = [int(x) for x in parts[2:6]]
                depth = sum(counts)
                
                if depth == 0:
                    continue
                
                # Get consensus base
                key = (contig, pos)
                if key not in consensus:
                    continue
                
                ref_base = consensus[key]['base']
                ref_idx = bases.index(ref_base) if ref_base in bases else -1
                
                # Find alternate alleles
                alt_bases = []
                for i, base in enumerate(bases):
                    if i != ref_idx:
                        if counts[i] >= args.min_obs and (counts[i] / depth) >= min_pct_decimal:
                            alt_bases.append(base)
                
                if not alt_bases:
                    continue
                
                alt_str = ','.join(alt_bases)
                
                # Check if SNP falls within a gene
                is_genic = False
                if contig in genes:
                    for gene in genes[contig]:
                        if gene['start'] <= pos <= gene['end']:
                            is_genic = True
                            # Calculate codon position
                            if gene['strand'] == 1:
                                pos_in_gene = pos - gene['start']
                                codon_num = pos_in_gene // 3 + 1
                                pos_in_codon = pos_in_gene % 3 + 1
                            else:
                                pos_in_gene = gene['end'] - pos
                                codon_num = pos_in_gene // 3 + 1
                                pos_in_codon = pos_in_gene % 3 + 1
                            
                            print(f"{contig}\t{pos}\t{ref_base}\t{alt_str}\t{depth}\t{gene['gene_id']}\t{gene['start']}\t{gene['end']}\t{gene['strand']}\t{codon_num}\t{pos_in_codon}\t{sample}", file=genic_out)
                
                if not is_genic:
                    print(f"{contig}\t{pos}\t{ref_base}\t{alt_str}\t{depth}\t{sample}", file=non_genic_out)
    
    genic_out.close()
    non_genic_out.close()
    
    # Write corrected genes (placeholder - would need full implementation)
    with open(args.output_genes, 'w') as f:
        with open(args.genes, 'r') as g:
            f.write(g.read())
    
    # Write codon counts (placeholder)
    with open(args.output_codons, 'w') as f:
        print("gene_id\t" + "\t".join([f"{b1}{b2}{b3}" for b1 in 'ATCG' for b2 in 'ATCG' for b3 in 'ATCG']), file=f)

if __name__ == "__main__":
    main()
