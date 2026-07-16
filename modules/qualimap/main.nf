#!/usr/bin/env nextflow

process QUALIMAP {
    label 'process_medium'

    input:
    tuple val(sample_name), path(bam)
    path(ref_gtf)

    output:
    path("${sample_name}_qualimap"), emit: qualimap_dir

    script:
    def mem = task.memory ? "${task.memory.toGiga()}G" : "16G"
    """
    qualimap rnaseq \\
        -bam ${bam} \\
        -gtf ${ref_gtf} \\
        -p strand-specific-reverse \\
        -outdir ${sample_name}_qualimap \\
        --java-mem-size="${mem}"
    """
}