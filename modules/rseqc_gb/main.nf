#!/usr/bin/env nextflow

process RSEQC_GB {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(ref_bed)

    output:
    path("*.txt"), emit: genebody
    path("*.pdf"), emit: plot, optional: true

    script:
    """
    geneBody_coverage.py -r ${ref_bed} -i ${bam} -o ${sample_name}
    """

}