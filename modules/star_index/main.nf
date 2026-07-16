#!/usr/bin/env nextflow

process STAR_INDEX {
    label 'process_high'
    clusterOptions '-l avx2'

    input:
    path fasta //params.reference_fasta
    path gtf //params.reference.gtf

    output:
    path("star_index"), emit: index // names the output as "index"

    script:
    """
    mkdir -p star

    STAR --runMode genomeGenerate \
        --genomeDir star_index \
        --genomeFastaFiles ${fasta} \
        --sjdbGTFfile ${gtf} \
        --runThreadN $task.cpus
    """
}