#!/usr/bin/env nextflow

// note: set -l flag to be average read length
// need to check the FASTQC report for those samples
// siRNA samples: ~51bp
// tissue samples: ~151bp
process RIBODETECTOR {
    label 'process_high'

    input:
    tuple val(sample_name), path(read1), path(read2)
    val(read_length)

    output:
    tuple val(sample_name), path("${sample_name}_norrna.1.fq"), path("${sample_name}_norrna.2.fq"), emit: filtered_fastqs
    path("*.log"), emit: logs

    script:
    """
    ribodetector_cpu -t ${task.cpus} \\
                    -l ${read_length} \\
                    -i ${read1} ${read2} \\
                    -e rrna \\
                    -o "${sample_name}_norrna.1.fq" "${sample_name}_norrna.2.fq" \\
                    > "${sample_name}_ribodetector.log"
    """
}