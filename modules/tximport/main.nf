#!/usr/bin/env nextflow

process TXIMPORT {
    label 'process_medium'

    input:
    path(quant_dirs)      // all salmon dirs, collected
    path(tx2gene)         // the .tsv, staged in

    output:
    path("gene_counts.tsv"), emit: gene_counts
    path("txi.rds"),         emit: txi

    script:
    """
    Rscript "${projectDir}/bin/tximport.R" ${tx2gene}
    """
}