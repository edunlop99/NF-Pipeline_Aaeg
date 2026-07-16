#!/usr/bin/env nextflow

process HTSEQ_COUNT {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam)
    path(gtf)

    output:
    tuple val(sample_name), path("${sample_name}.htseq.counts"), emit: counts

    script:
    """
    htseq-count \\
        -f bam \\
        -r pos \\
        -s reverse \\
        -t exon \\
        -i gene_id \\
        ${bam} \\
        ${gtf} > ${sample_name}.htseq.counts
    """
}