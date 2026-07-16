#!/usr/bin/env nextflow

process DEDUPLICATE {
    label 'process_high'

    input:
    tuple val(sample_name), path(bam_file)   // from STAR_ALIGN (coordinate-sorted BAM)

    output:
    tuple val(sample_name), path("${sample_name}.dedup.bam"), path("${sample_name}.dedup.bam.bai"), emit: dedup_bam
    tuple val(sample_name), path("${sample_name}_dedup.log"), emit: log

    script:
    """
    sambamba markdup -r -t ${task.cpus} \\
        ${bam_file} "${sample_name}.dedup.bam" \\
        2> "${sample_name}_dedup.log"

    samtools index -@ ${task.cpus} "${sample_name}.dedup.bam"
    """
}