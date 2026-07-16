#!/usr/bin/env nextflow

// Using Zymo-Seq Tibo-Free Total RNA adapter sequences:
// -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC 
// -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGA

process TRIM_GALORE {
    label 'process_medium'

    input:
    tuple val(sample_name), path(read1), path(read2)

    output:
    tuple val(sample_name), path("${sample_name}_val_1.fq*"), path("${sample_name}_val_2.fq*"), emit: trimmed_fastqs
    path("*_trimming_report.txt"), emit: reports

    script:
    """
    trim_galore --cores ${task.cpus} \\
                --paired --basename ${sample_name} \\
                -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \\
                -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGA \\
                ${read1} ${read2}
    """
}