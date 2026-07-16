#!/usr/bin/env nextflow

process BOWTIE2_INDEX {
    label 'process_medium'

    input:
    path rrna_fasta
    path virna_fasta

    output:
    path "rrna_index", emit: rrna_index
    path "virna_index", emit: virna_index

    script:
    """
    mkdir rrna_index virna_index
    bowtie2-build --threads ${task.cpus} ${rrna_fasta} rrna_index/rrna
    bowtie2-build --threads ${task.cpus} ${virna_fasta} virna_index/virna
    """
}