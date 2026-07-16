#!/usr/bin/env nextflow

process DESEQ2 {
    label 'process_low'

    input:
    path(featurecounts_matrix)      // featureCounts raw counts output
    path(salmon_txi)                // Salmon quant txi.rds

    output:
    path("deseq2_outdir"), emit: outdir

    script:
    """
    Rscript "${projectDir}/bin/deseq2.R" ${featurecounts_matrix} ${salmon_txi} "deseq2_outdir"
    """
}