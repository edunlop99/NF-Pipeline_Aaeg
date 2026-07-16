#!/usr/bin/env nextflow

process SALMON_QUANT {
    label 'process_high'

    input:
    tuple val(sample_name), path(read1), path(read2)
    path(salmon_index)

    output:
    tuple val(sample_name), path("${sample_name}_salmon"), emit: quant_dir

    script:
    """
    salmon quant -i ${salmon_index} \\
                -l ISR \\
                -1 ${read1} -2 ${read2} \\
                -p ${task.cpus} \\
                --validateMappings \\
                -o "${sample_name}_salmon"
    """

}