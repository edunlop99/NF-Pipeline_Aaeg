#!/usr/bin/env nextflow

process STAR_ALIGN {
    label 'process_high'
    clusterOptions '-l avx2'

    input:
    tuple val(sample_name), path(reads_fwd), path(reads_rev)
    path(index)

    output:
    tuple val(sample_name), path("${sample_name}_Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(sample_name), path("${sample_name}_ReadsPerGene.out.tab"), emit: gene_counts
    tuple val(sample_name), path("${sample_name}_Log.final.out"), emit: log

    script: 
    def gz = reads_fwd.name.endsWith('.gz') ? '--readFilesCommand zcat' : ''
    """ 
    STAR --runThreadN ${task.cpus} \\
        --genomeDir ${index} \\
        --readFilesIn ${reads_fwd} ${reads_rev} \\
        ${gz} \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --outFileNamePrefix "./${sample_name}_"
    """
}
