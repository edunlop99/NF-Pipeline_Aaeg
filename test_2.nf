#!/usr/bin/env nextflow

// Define Processes
include {FASTQC} from './modules/fastqc'
include {MULTIQC} from './modules/multiqc'
include {TRIM_GALORE} from './modules/trim_galore'
include {RIBODETECTOR} from './modules/ribodetector'
include {VIRAL_COUNTS} from './modules/viral_counts'
include {VIRAL_MATRIX} from './modules/viral_matrix'
include {BOWTIE2_INDEX} from './modules/bowtie2_index'
include {BOWTIE2_ALIGN} from './modules/bowtie2_align'
include {STAR_INDEX} from './modules/star_index'
include {STAR_ALIGN} from './modules/star_align'
include {SALMON_INDEX} from './modules/salmon_index'
include {SALMON_QUANT} from './modules/salmon_quant'
include {TXIMPORT} from './modules/tximport'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {SAMTOOLS_INDEX} from './modules/samtools_index'
include {DEDUPLICATE} from './modules/deduplicate'
include {FEATURECOUNTS} from './modules/feature_counts'
include {HTSEQ_COUNT} from './modules/htseq_count'
include {STRINGTIE} from './modules/stringtie'
include {PREPDE} from './modules/prep_de'
include {BAMCOVERAGE} from './modules/deeptools_bamcoverage'
include {QUALIMAP} from './modules/qualimap'

workflow QC_only {
    // just for running QC on raw FastQ files
    take:
    read_pairs_ch

    main:
    fastqc_ch = read_pairs_ch.flatMap { sample_name, fastq_read1, fastq_read2 ->  
        [[sample_name, fastq_read1], [sample_name, fastq_read2]]  
    }

    // FastQC and MultiQC for raw FASTQs
    FASTQC(fastqc_ch)
    fastqc_out_html = FASTQC.out.html
    fastqc_zip_ch = FASTQC.out.zip
    fastqc_zip_ch.map{ it[1] }.collect()
    | set { multiqc_ch }
    multiqc_report_01 = MULTIQC( multiqc_ch )

    emit: 
    multiqc_report_01

}

workflow STAR_alignment {
    take:
    read_pairs_ch

    main:
    // Adapter Trimming
    TRIM_GALORE(read_pairs_ch)
    trimmed_files_ch = TRIM_GALORE.out.trimmed_fastqs
    trim_reports = TRIM_GALORE.out.reports

    // rRNA alignment and viral RNA alignment
    BOWTIE2_INDEX(rrna_fasta, virna_fasta)
    BOWTIE2_ALIGN(trimmed_files_ch, BOWTIE2_INDEX.out.rrna_index, BOWTIE2_INDEX.out.virna_index)
    rrna_virna_logs = BOWTIE2_ALIGN.out.logs

    // Get viral RNA counts
    VIRAL_COUNTS(BOWTIE2_ALIGN.out.virna_bam) // tuple val(sample_name), path(bam) 
    VIRAL_COUNTS.out.counts.map { name, tsv -> tsv }.collect()
    | set { viral_counts_ch }
    viral_counts_matrix = VIRAL_MATRIX(viral_counts_ch)
    
    // Ribodetector; ML-based ribosomal RNA detection, output filtered fastq's
    RIBODETECTOR(trimmed_files_ch, read_length)
    filtered_fastqs = RIBODETECTOR.out.filtered_fastqs
    ribodetector_logs = RIBODETECTOR.out.logs

    // Generate STAR Index and align
    STAR_INDEX(ref_genome, ref_gtf)
    STAR_ALIGN(filtered_fastqs, STAR_INDEX.out.index) // output bam files tuple val(sample_name), path(bam)
    star_align_out = STAR_ALIGN.out.bam
    STAR_ALIGN.out.log.map{ it[1] }.collect() // collect all STAR output logs into list
    | set { star_logs } 
    STAR_ALIGN.out.gene_counts.map{ it[1] }.collect()
    | set { star_gene_counts }
    
    // Post-alignment MultiQC
    trim_reports
        .mix(star_logs, star_gene_counts, rrna_virna_logs, ribodetector_logs)
        .collect()
        | set { post_align_multiqc_ch }
    post_align_multiqc = MULTIQC(post_align_multiqc_ch)

    emit:
    rrna_virna_logs = rrna_virna_logs
    viral_counts_matrix = viral_counts_matrix
    filtered_fastqs = filtered_fastqs
    star_align_out = star_align_out
    post_align_multiqc = post_align_multiqc

}

workflow Post_alignment_processing {
    take:
    star_align_out

    main:
    // Index and Deduplication
    // NOTE: STAR --outSAMtype SoortedByCoordinate means that STAR_ALIGN.out.bam is already sorted, SAMTOOLS sort is redundant
    DEDUPLICATE(star_align_out) // note: doesn't actually remove duplicates, just labels them as duplicates
    dedup_ch = DEDUPLICATE.out.dedup_bam   // (sample, dedup.bam, dedup.bam.bai)
    dedup_logs = DEDUPLICATE.out.log.map { it[1] }

    // Coverage Plots
    // input for bam coverage: tuple val(sample_name), path(bam), path(bai)
    bigwig_out= BAMCOVERAGE(dedup_ch)

    // Quality control
    dedup_bams_ch = dedup_ch.map { it -> tuple(it[0], it[1]) }
    QUALIMAP(dedup_bams_ch, ref_gtf)
    qualimap_stats = QUALIMAP.out.qualimap_dir
    SAMTOOLS_FLAGSTAT(dedup_bams_ch)
    samtools_stats = SAMTOOLS_FLAGSTAT.out.flagstat

    // Counts Extraction
    all_bams_ch = dedup_ch.map { it -> it[1] }.collect()
    FEATURECOUNTS(all_bams_ch, ref_gtf)
    feature_counts = FEATURECOUNTS.out.counts
    feature_counts_summary = FEATURECOUNTS.out.summary
    //HTSEQ_COUNT(dedup_bams_ch, ref_gtf)
    //htseq_counts = HTSEQ_COUNT.out.counts // tuple val(sample_name), path(counts)

    SALMON_INDEX(ref_transcriptome)
    SALMON_QUANT(dedup_bams_ch, SALMON_INDEX.out.index) // transcriptome-based quantification; note that this will have discrepency with FeatureCounts counts
    SALMON_QUANT.out.quant_dir
        .map { name, dir -> dir }
        .collect()
        | set { salmon_dirs }

    TXIMPORT(salmon_dirs)
    salmon_counts = TXIMPORT.out.gene_counts

    // Stringtie Transcript Quantification
    // NOTE: StringTie fails due to lines in NCBI GTF with 'transcript_id "";' so fix is to remove these lines, which also means all the pipeline processes use this GTF instead.
    // STRINGTIE(DEDUPLICATE.out.dedup_bam, ref_gtf)
    // STRINGTIE.out.quant_gtf
    //     .multiMap { name, gtf ->
    //         names: name
    //         gtfs:  gtf
    //     }
    //     .set { collected }

    // PREPDE(collected.gtfs.collect(), collected.names.collect())
    // stringtie_counts = PREPDE.out.gene_counts

    // MultiQC
    qualimap_stats
        .mix(samtools_stats, dedup_logs, feature_counts_summary, salmon_dirs)
        .collect()
        | set { final_multiqc_ch }
    final_multiqc_report = MULTIQC(final_multiqc_ch)

    emit:
    dedup_ch // (sample, dedup.bam, dedup.bam.bai)
    bigwig_out
    feature_counts
    salmon_counts
    //stringtie_counts
    final_multiqc_report
}

workflow {

    main:
    // define params inputs
    reads = params.reads
    read_length = params.read_length
    ref_genome = params.reference_genome
    ref_gtf = params.reference_gtf
    ref_transcriptome = params.reference_transcriptome
    rrna_fasta = params.rrna_fasta
    virna_fasta = params.virna_fasta

    // Parse sample sheet
    Channel.fromPath(params.reads)
    | splitCsv(header:true)
    | map { row -> tuple(row.sample_name, row.fastq_read1, row.fastq_read2) }
    | set { read_pairs_ch } // tuple val(name), path(read 1), path(read 2)

    // FastQC/MultiQC for raw FastQ's
    QC_only( read_pairs_ch )

    // Trimming, rRNA/viRNA count, STAR, alignment QC
    STAR_alignment( read_pairs_ch )

    // Generate counts, BAM Coverage plots, final QC
    Post_alignment_processing()

    publish:
    fastqc_out_html = fastqc_out_html
    fastqc_out_zip = fastqc_zip_ch
    multiqc_report_01 = multiqc_report_01
    trim_galore_out = trimmed_files_ch
    rrna_virna_logs = rrna_virna_logs
    star_align_out = star_align_out
    post_align_multiqc = post_align_multiqc
    viral_counts_matrix = viral_counts_matrix
    bigwig_out = bigwig_out
    feature_counts = feature_counts
    //htseq_counts = htseq_counts
    //stringtie_counts = stringtie_counts_ch
    salmon_counts = salmon_counts
    final_multiqc_report = final_multiqc_report
}