#!/usr/bin/env nextflow

process SAMTOOLS_FLAGSTAT {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam)

    output:
    path("*flagstat.txt"), emit: flagstat

    shell:
    """
    samtools flagstat -@ $task.cpus $bam > "${sample_name}_flagstat.txt"
    """
}