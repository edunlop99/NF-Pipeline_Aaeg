#!/usr/bin/env nextflow

process STRINGTIE {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(ref_gtf)

    output:
    tuple val(sample_name), path("${sample_name}/${sample_name}.gtf"), emit: quant_gtf
    tuple val(sample_name), path("${sample_name}_gene_abund.tab"), emit: gene_abundance

    script:
    """
    stringtie ${bam} \\
        -e \\
        -G ${ref_gtf} \\
        -o ${sample_name}/${sample_name}.gtf \\
        -A ${sample_name}_gene_abund.tab \\
        -p ${task.cpus} \\
        --rf
    """
}