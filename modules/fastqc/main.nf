#!/usr/bin/env nextflow

process FASTQC {
    label 'process_low'

    input:
    tuple val(sample_name), path(fastq)  // Input value is a tuple, including sample name and the paths to the files

    output:
    tuple val(sample_name), path("*.zip"), emit: zip // first output, labeled "zip"
    tuple val(sample_name), path('*.html'), emit: html // second output, labeled "html"

    script:
    """
    fastqc -t $task.cpus $fastq
    """
}