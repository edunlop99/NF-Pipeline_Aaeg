#!/usr/bin/env nextflow

process RSEQC_RD {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(ref_bed)

    output:
    tuple val(sample_name), path("${sample_name}_read_distribution.txt"), emit: readdist

    script:
    """
    read_distribution.py  -i ${bam} -r ${params.reference_bed} > "${sample_name}_read_distribution.txt"

    """
}