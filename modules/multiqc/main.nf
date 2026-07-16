#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_low'

    input:
    path('*') //anything in current directory

    output:
    path("*.html") // expect .html

    script:
    """
    multiqc .
    """
}