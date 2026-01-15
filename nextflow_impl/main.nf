#!/usr/bin/env nextflow

// =============================================================================
// METAPOP - Nextflow Pipeline for Metagenomic Population Analysis
// =============================================================================
// A pipeline for macro- and micro-diversity analyses of metagenomic-derived 
// populations. Reimplemented from the original MetaPop Python/R pipeline.
// =============================================================================

nextflow.enable.dsl = 2

// =============================================================================
// Import Modules
// =============================================================================

include { 
    INDEX_REFERENCE;
    PREDICT_GENES;
    FILL_MD; 
    FILTER_READS; 
    CALCULATE_COVERAGE; 
    FILTER_GENOMES 
} from './modules/preprocessing'

include { 
    CALL_VARIANTS;
    REFINE_SNPS;
    CALL_CONSENSUS;
    ASSOCIATE_GENES
} from './modules/snp_calling'

include {
    MINE_READS;
    CALCULATE_MICRODIVERSITY;
    CALCULATE_FST
} from './modules/microdiversity'

// =============================================================================
// Validate Inputs
// =============================================================================

if (params.input_dir == null) {
    error "Please provide input BAM directory with --input_dir"
}

if (params.reference == null) {
    error "Please provide reference genome with --reference"
}

// =============================================================================
// Main Workflow
// =============================================================================

workflow {
    
    log.info """
    ╔═══════════════════════════════════════════════════════════════════════════╗
    ║                                                                           ║
    ║   M E T A P O P   -   N E X T F L O W   P I P E L I N E                  ║
    ║                                                                           ║
    ╠═══════════════════════════════════════════════════════════════════════════╣
    ║   Input BAMs    : ${params.input_dir}
    ║   Reference     : ${params.reference}
    ║   Genes         : ${params.genes ?: 'Will be predicted'}
    ║   Output        : ${params.output_dir}
    ╚═══════════════════════════════════════════════════════════════════════════╝
    """

    // -------------------------------------------------------------------------
    // Setup Channels
    // -------------------------------------------------------------------------
    
    // Input BAM files
    input_bams_ch = Channel.fromPath("${params.input_dir}/*.bam")
        .map { bam -> 
            def meta = [id: bam.simpleName]
            tuple(meta, bam)
        }
    
    // Reference genome
    reference_ch = Channel.fromPath(params.reference)
    
    // -------------------------------------------------------------------------
    // STAGE 0: Reference Preparation
    // -------------------------------------------------------------------------
    
    // Index reference
    INDEX_REFERENCE(reference_ch)
    
    // Predict genes (if not provided)
    if (params.genes) {
        genes_ch = Channel.fromPath(params.genes)
    } else {
        PREDICT_GENES(INDEX_REFERENCE.out.fasta)
        genes_ch = PREDICT_GENES.out.genes
    }
    
    // -------------------------------------------------------------------------
    // STAGE 1: Preprocessing
    // -------------------------------------------------------------------------
    
    if (!params.skip_preprocessing) {
        // Fill MD tags for percent identity calculation
        FILL_MD(input_bams_ch, INDEX_REFERENCE.out.fasta)
        
        // Filter reads by length and percent identity
        FILTER_READS(FILL_MD.out.bam)
        
        // Calculate breadth and depth of coverage
        CALCULATE_COVERAGE(FILTER_READS.out.bam, INDEX_REFERENCE.out.fai)
        
        // Join filtered BAMs with passing genomes list
        filter_input_ch = FILTER_READS.out.bam
            .join(CALCULATE_COVERAGE.out.pass_list)
        
        // Filter to keep only passing genomes
        FILTER_GENOMES(filter_input_ch)
        
        preprocessed_bams_ch = FILTER_GENOMES.out.bam
    } else {
        // If skipping preprocessing, use input BAMs directly
        preprocessed_bams_ch = input_bams_ch
            .map { meta, bam -> tuple(meta, bam, file("${bam}.bai")) }
    }
    
    // -------------------------------------------------------------------------
    // STAGE 2: SNP Calling
    // -------------------------------------------------------------------------
    
    if (!params.skip_snp_calling) {
        // Call variants with bcftools
        CALL_VARIANTS(preprocessed_bams_ch, INDEX_REFERENCE.out.fasta)
        
        // Combine all variant positions
        all_variants_ch = CALL_VARIANTS.out.variants
            .map { meta, variants -> variants }
            .collectFile(name: 'all_variants.txt', newLine: true)
        
        // Refine SNP calls
        REFINE_SNPS(
            preprocessed_bams_ch, 
            INDEX_REFERENCE.out.fasta, 
            all_variants_ch
        )
        
        // Collect all SNP count files for consensus calling
        snp_counts_ch = REFINE_SNPS.out.snp_counts
            .map { meta, counts -> counts }
            .collect()
        
        // Call consensus bases
        CALL_CONSENSUS(snp_counts_ch)
        
        // Associate SNPs with genes
        ASSOCIATE_GENES(
            snp_counts_ch,
            CALL_CONSENSUS.out.consensus,
            genes_ch,
            INDEX_REFERENCE.out.fasta
        )
    }
    
    // -------------------------------------------------------------------------
    // STAGE 3: Microdiversity (if implemented)
    // -------------------------------------------------------------------------
    
    // Placeholder for microdiversity analysis
    // This would include:
    // - MINE_READS for linked SNP analysis
    // - CALCULATE_MICRODIVERSITY for Pi, Theta, pN/pS
    // - CALCULATE_FST for population differentiation
    
    // -------------------------------------------------------------------------
    // Summary
    // -------------------------------------------------------------------------
}

// =============================================================================
// Workflow Completion Handler
// =============================================================================

workflow.onComplete {
    log.info """
    ============================================================
                    Pipeline Complete!
    ============================================================
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    Results     : ${params.output_dir}
    ============================================================
    """
}
