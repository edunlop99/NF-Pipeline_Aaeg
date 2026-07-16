#!/usr/bin/env nextflow

process PLOT_DEGS {
    label 'process_low'

    input:
    path(fc_deseq2_out) // featureCounts deseq2 dds.rds
    path(slmn_deseq2_out) // salmon counts deseq2 dds.rds

    output:
    path("*.png"), emit: plots

    script:
    """
    Rscript "${projectDir}/bin/plot_degs.R" ${fc_deseq2_out} ${slmn_deseq2_out}
    """
}