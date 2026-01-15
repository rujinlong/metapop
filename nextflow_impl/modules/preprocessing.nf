// =============================================================================
// PREPROCESSING MODULE
// =============================================================================
// Implements the preprocessing steps of MetaPop:
// 1. MD tagging for percent identity calculation
// 2. Read filtering by length and percent identity
// 3. Coverage calculation (breadth and TAD)
// 4. Genome filtering based on coverage thresholds
// =============================================================================

process FILL_MD {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.output_dir}/preprocessing/md_tagged", mode: 'symlink'

    input:
    tuple val(meta), path(bam)
    path reference

    output:
    tuple val(meta), path("${meta.id}_md.bam"), path("${meta.id}_md.bam.bai"), emit: bam

    script:
    """
    samtools calmd -b $bam $reference > ${meta.id}_md.bam
    samtools index ${meta.id}_md.bam
    """
}

process FILTER_READS {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.output_dir}/preprocessing/filtered_reads", mode: 'symlink'

    input:
    tuple val(meta), path(bam), path(bai)
    
    output:
    tuple val(meta), path("${meta.id}_filtered.bam"), path("${meta.id}_filtered.bam.bai"), emit: bam

    script:
    def global_flag = params.is_global ? "--global_id" : ""
    """
    filter_reads.py \\
        --input $bam \\
        --output ${meta.id}_filtered.bam \\
        --min_len ${params.min_len} \\
        --min_pct_id ${params.min_pct_id} \\
        $global_flag
        
    samtools index ${meta.id}_filtered.bam
    """
}

process CALCULATE_COVERAGE {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.output_dir}/preprocessing/coverage", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta_index

    output:
    tuple val(meta), path("${meta.id}_breadth_depth.tsv"), emit: stats
    tuple val(meta), path("${meta.id}_passing_genomes.txt"), emit: pass_list
    tuple val(meta), path("${meta.id}_depth_per_pos.tsv.gz"), emit: depth_per_pos

    script:
    """
    # Create lengths file from fai
    cut -f1,2 $fasta_index > lengths.tsv

    # Calculate depth per position and pipe to coverage calculator
    samtools depth -a $bam | tee >(gzip > ${meta.id}_depth_per_pos.tsv.gz) | \\
        calc_coverage.py \\
            --lengths lengths.tsv \\
            --min_cov ${params.min_cov} \\
            --min_dep ${params.min_dep} \\
            --trunc ${params.truncation} \\
            --output_stats ${meta.id}_breadth_depth.tsv \\
            --output_pass ${meta.id}_passing_genomes.txt
    """
}

process FILTER_GENOMES {
    tag "$meta.id"
    label 'process_medium'
    
    publishDir "${params.output_dir}/preprocessing/preprocessed_bams", mode: 'symlink'

    input:
    tuple val(meta), path(bam), path(bai), path(pass_list)

    output:
    tuple val(meta), path("${meta.id}_preprocessed.bam"), path("${meta.id}_preprocessed.bam.bai"), emit: bam

    script:
    """
    filter_genomes.py \\
        --input $bam \\
        --output ${meta.id}_preprocessed.bam \\
        --genomes $pass_list
        
    samtools index ${meta.id}_preprocessed.bam
    """
}

process INDEX_REFERENCE {
    label 'process_low'
    
    publishDir "${params.output_dir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "${fasta}.fai", emit: fai
    path fasta, emit: fasta

    script:
    """
    samtools faidx $fasta
    """
}

process PREDICT_GENES {
    label 'process_medium'
    
    // If single mode fails (sequences too small), retry with meta mode
    errorStrategy { task.attempt <= 1 && params.pyrodigal_mode == 'single' ? 'retry' : 'terminate' }
    maxRetries 1
    
    publishDir "${params.output_dir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "genes.fna", emit: genes
    path "genes.faa", emit: proteins
    path "genes.gff", emit: gff

    script:
    // Use meta mode on retry (for small sequences that can't be trained)
    def mode = task.attempt > 1 ? 'meta' : params.pyrodigal_mode
    """
    pyrodigal -i $fasta \\
        -a genes.faa \\
        -d genes.fna \\
        -o genes.gff \\
        -p ${mode} \\
        -f gff \\
        -m \\
        --pool process \\
        -j ${task.cpus}
    """
}
