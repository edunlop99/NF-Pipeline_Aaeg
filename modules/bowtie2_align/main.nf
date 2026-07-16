#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    // Process to align fastq's to viral RNA reference w/ Bowtie2
    label 'process_high'
    errorStrategy 'finish'

    input:
    tuple val(sample_name), path(reads_fwd), path(reads_rev)
    path rrna_index
    path virna_index

    output:
    tuple val(sample_name), path("${sample_name}_rrna.sorted.bam"), emit: rrna_bam 
    tuple val(sample_name), path("${sample_name}_virna.sorted.bam"), emit: virna_bam 
    path('*log'), emit: logs

    script:
    """
    bowtie2 -x rrna_index/rrna \\
        -1 ${reads_fwd} -2 ${reads_rev} \\
        -p ${task.cpus} \\
        2> ${sample_name}_rrna_align.log \\
        | samtools view -bS -F 4 \\
        | samtools sort -o ${sample_name}_rrna.sorted.bam

    bowtie2 -x virna_index/virna \\
        -1 ${reads_fwd} -2 ${reads_rev} \\
        -p ${task.cpus} \\
        2> ${sample_name}_virna_align.log \\
        | samtools view -bS -F 4 \\
        | samtools sort -o ${sample_name}_virna.sorted.bam
    """

}