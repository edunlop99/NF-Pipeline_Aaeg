#!/usr/bin/env nextflow

process RSEQC_GB {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam), path(bai)
    path(ref_bed)

    output:
    path("${sample_name}.geneBodyCoverage.txt"), emit: genebody
    path("${sample_name}.geneBodyCoverage.curves.pdf"), emit: plot, optional: true

    script:
    """
    geneBody_coverage.py -r ${ref_bed} -i ${bam} -o ${sample_name}
    """

}