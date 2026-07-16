#!/usr/bin/env nextflow

process SALMON_INDEX {
    label 'process_high'

    input:
    path(transcriptome_fasta)

    output:
    path "salmon_index", emit: index

    script:
    """
    salmon index -t ${transcriptome_fasta} -i salmon_index -p ${task.cpus}
    """

}