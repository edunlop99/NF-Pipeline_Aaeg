#!/usr/bin/env nextflow

process PREPDE {
    label 'process_low'

    input:
    path(gtf_files)          // all per-sample GTFs, collected
    val(sample_names)        // matching sample names

    output:
    path("gene_count_matrix.csv"), emit: gene_counts
    path("transcript_count_matrix.csv"), emit: transcript_counts

    script:
    """
    # build the sample list file prepDE.py expects: <sample_name> <path_to_gtf>
    paste -d' ' \\
        <(printf '%s\\n' ${sample_names.join(' ')}) \\
        <(ls -1 *.gtf) > sample_list.txt

    prepDE.py -i sample_list.txt \\
        -g gene_count_matrix.csv \\
        -t transcript_count_matrix.csv
    """
}