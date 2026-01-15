// =============================================================================
// SNP CALLING MODULE
// =============================================================================
// Implements the SNP calling steps of MetaPop:
// 1. Variant calling with bcftools
// 2. SNP refinement with samtools mpileup
// 3. Consensus base calling
// 4. Gene association
// =============================================================================

process CALL_VARIANTS {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.output_dir}/snp_calling/variants", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${meta.id}_variants.txt"), emit: variants

    script:
    """
    # Create ploidy file for haploid organisms
    echo -e "*\\t*\\t*\\t*\\t1" > ploidy.txt
    
    # Call variants with bcftools
    bcftools mpileup -Ob -I -f $reference -o ${meta.id}.mpileup.bcf $bam
    bcftools call -mv -Ob --ploidy-file ploidy.txt -o ${meta.id}.calls.bcf ${meta.id}.mpileup.bcf
    bcftools filter -i "QUAL>${params.min_qual}" -Ob -o ${meta.id}.filtered.bcf ${meta.id}.calls.bcf
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' -o ${meta.id}_variants.txt ${meta.id}.filtered.bcf
    
    # Cleanup intermediate files
    rm -f ${meta.id}.mpileup.bcf ${meta.id}.calls.bcf ${meta.id}.filtered.bcf
    """
}

process REFINE_SNPS {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.output_dir}/snp_calling/snp_counts", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path all_variants

    output:
    tuple val(meta), path("${meta.id}_snp_counts.tsv"), emit: snp_counts

    script:
    """
    refine_snps.py \\
        --bam $bam \\
        --reference $reference \\
        --variants $all_variants \\
        --min_qual ${params.min_qual} \\
        --output ${meta.id}_snp_counts.tsv
    """
}

process CALL_CONSENSUS {
    label 'process_low'
    
    publishDir "${params.output_dir}/snp_calling/consensus", mode: 'copy'

    input:
    path snp_counts  // All SNP count files

    output:
    path "consensus_bases.tsv", emit: consensus

    script:
    """
    call_consensus.py \\
        --input_dir . \\
        --output consensus_bases.tsv
    """
}

process ASSOCIATE_GENES {
    label 'process_medium'
    
    publishDir "${params.output_dir}/snp_calling/genic_snps", mode: 'copy'

    input:
    path snp_counts  // All SNP count files
    path consensus
    path genes
    path reference

    output:
    path "genic_snps.tsv", emit: genic_snps
    path "non_genic_snps.tsv", emit: non_genic_snps
    path "corrected_genes.fna", emit: corrected_genes
    path "codon_counts.tsv", emit: codon_counts

    script:
    """
    associate_genes.py \\
        --snp_dir . \\
        --consensus $consensus \\
        --genes $genes \\
        --reference $reference \\
        --min_obs ${params.min_obs} \\
        --min_pct ${params.min_pct} \\
        --output_genic genic_snps.tsv \\
        --output_non_genic non_genic_snps.tsv \\
        --output_genes corrected_genes.fna \\
        --output_codons codon_counts.tsv
    """
}
