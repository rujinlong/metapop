// =============================================================================
// MICRODIVERSITY MODULE
// =============================================================================
// Implements the microdiversity analysis steps of MetaPop:
// 1. SNP linking (linked SNPs on same codon)
// 2. Microdiversity statistics (Pi, Theta, Tajima's D, pN/pS)
// 3. FST calculation
// =============================================================================

process MINE_READS {
    tag "all_samples"
    label 'process_high'
    
    publishDir "${params.output_dir}/microdiversity/linked_snps", mode: 'copy'

    input:
    path genic_snps
    path preprocessed_bams  // All BAM files

    output:
    path "linked_snp_results.tsv", emit: linked_snps

    script:
    """
    mine_reads.py \\
        --snps $genic_snps \\
        --bam_dir . \\
        --output linked_snp_results.tsv
    """
}

process CALCULATE_MICRODIVERSITY {
    label 'process_medium'
    
    publishDir "${params.output_dir}/microdiversity", mode: 'copy'

    input:
    path genic_snps
    path linked_snps
    path genes
    path reference
    path fai

    output:
    path "global_microdiversity.tsv", emit: global_stats
    path "per_gene_microdiversity.tsv", emit: per_gene_stats
    path "selection_pressure.tsv", emit: selection

    script:
    """
    calculate_microdiversity.py \\
        --genic_snps $genic_snps \\
        --linked_snps $linked_snps \\
        --genes $genes \\
        --reference $reference \\
        --fai $fai \\
        --subsample ${params.subsample_size} \\
        --output_global global_microdiversity.tsv \\
        --output_genes per_gene_microdiversity.tsv \\
        --output_selection selection_pressure.tsv
    """
}

process CALCULATE_FST {
    label 'process_medium'
    
    publishDir "${params.output_dir}/microdiversity", mode: 'copy'

    input:
    path microdiversity_data
    path contig_lengths

    output:
    path "fst_results.tsv", emit: fst

    script:
    """
    calculate_fst.py \\
        --microdiversity $microdiversity_data \\
        --lengths $contig_lengths \\
        --output fst_results.tsv
    """
}
