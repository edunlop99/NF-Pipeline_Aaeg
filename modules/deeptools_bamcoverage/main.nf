#!/usr/bin/env nextflow

process BAMCOVERAGE {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam), path(bai)

    output:
    tuple val(sample_name), path("${sample_name}_fwd.bw"), path("${sample_name}_rev.bw"), emit: bigwigs

    shell:
    """
    bamCoverage -b ${bam} -o "${sample_name}_fwd.bw" --filterRNAstrand forward -p $task.cpus
    bamCoverage -b ${bam} -o "${sample_name}_rev.bw" --filterRNAstrand reverse -p $task.cpus
    """
}
