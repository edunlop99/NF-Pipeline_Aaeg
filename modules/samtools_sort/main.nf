#!/usr/bin/env nextflow

process SAMTOOLS_SORT {
    label 'process_high'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam"), emit: bam_sorted

    script:
    """ 
    samtools sort -@ ${task.cpus} -o "${sample_name}.sorted.bam" ${bam}
    """

}