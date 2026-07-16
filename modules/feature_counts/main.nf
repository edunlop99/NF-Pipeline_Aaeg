#!/usr/bin/env nextflow

process FEATURECOUNTS {
    label 'process_medium'

    input:
    path(bams)
    path(gtf)

    output:
    path("counts_matrix.txt"), emit: counts
    path("counts_matrix.txt.summary"), emit: summary

    script:
    """
    featureCounts -T ${task.cpus} -p -M \\
        --countReadPairs -s 2 \\
        -t exon -g gene_id -a ${gtf} \\
        -o counts_matrix.txt ${bams}

    """
}