#!/usr/bin/env nextflow

process SAMTOOLS_INDEX {
    label 'process_medium'

    input:
    tuple val(sample_name), path(sorted_bam) // tuple of name, sorted BAM file

    output:
    tuple val(sample_name), path(sorted_bam), path("${sample_name}.bai"), emit: bai

    script:
    """
    samtools index --@ ${task.cpus} -o "${sample_name}.bai"  ${sorted_bam}
    """

}