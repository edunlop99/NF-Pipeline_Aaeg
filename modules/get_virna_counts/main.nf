#!/usr/bin/env nextflow

process GET_VIRNA_COUNTS {
    label 'process_medium'
    errorStrategy 'finish'

    input:
    tuple val(sample_name), path(virna_bam) // viral RNA alignment BAM (Bowtie2 output)

    output:
    path(virna_raw_counts)

    script:
    """
    samtools index ${virna_bam} \\
    samtools idxstats ${virna_bam}
    
    """
}