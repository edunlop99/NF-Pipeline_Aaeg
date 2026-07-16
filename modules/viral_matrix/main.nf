#!/usr/bin/env nextflow

process VIRAL_MATRIX {
    label 'process_single'

    input:
    path(count_files)

    output:
    path("viral_counts_matrix.tsv")

    script:
    """
    # build header: 'virus' + one column per sample
    header="virus"
    for f in *_viral_counts.tsv; do
        s=\$(basename "\$f" _viral_counts.tsv)
        header="\${header}\\t\${s}"
    done
    echo -e "\$header" > viral_counts_matrix.tsv

    # reference names from first file
    first=\$(ls *_viral_counts.tsv | head -1)
    cut -f1 "\$first" > refs.txt

    # paste all count columns beside the names
    cols=refs.txt
    for f in *_viral_counts.tsv; do
        cut -f2 "\$f" > "\${f}.col"
        cols="\$cols \${f}.col"
    done
    paste \$cols >> viral_counts_matrix.tsv
    """
}