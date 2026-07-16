#!/usr/bin/env nextflow

process VIRAL_COUNTS {
    label 'process_low'

    input:
    tuple val(sample_name), path(bam)

    output:
    tuple val(sample_name), path("${sample_name}_viral_counts.tsv"), emit: counts

    script:
    """
    samtools index ${bam}
    samtools idxstats ${bam} > ${sample_name}_viral_idxstats.tsv

    # idxstats columns: refname, seqlength, mapped, unmapped
    # keep the viral refs (drop the '*' unmapped line), emit refname + mapped count
    grep -v '^\\*' ${sample_name}_viral_idxstats.tsv | cut -f1,3 > ${sample_name}_viral_counts.tsv
    """
}