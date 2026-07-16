#!/usr/bin/env nextflow

// Define Processes
include {FASTQC} from './modules/fastqc'
include {MULTIQC as MULTIQC_01} from './modules/multiqc'
include {TRIM_GALORE} from './modules/trim_galore'
include {RIBODETECTOR} from './modules/ribodetector'
include {VIRAL_COUNTS} from './modules/viral_counts'
include {VIRAL_MATRIX} from './modules/viral_matrix'
include {BOWTIE2_INDEX} from './modules/bowtie2_index'
include {BOWTIE2_ALIGN} from './modules/bowtie2_align'
include {STAR_INDEX} from './modules/star_index'
include {STAR_ALIGN} from './modules/star_align'
include {MULTIQC as MULTIQC_02} from './modules/multiqc'
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
include {RSEQC_RD} from './modules/rseqc_rd'
include {RSEQC_GB} from './modules/rseqc_gb'
include {MULTIQC as MULTIQC_03} from './modules/multiqc'
include {DESEQ2} from './modules/deseq2'

workflow  {

    main:
    // define inputs
    reads = params.reads
    ref_genome = params.reference_genome
    ref_gtf = params.reference_gtf
    ref_transcriptome = params.reference_transcriptome
    ref_bed = params.reference_bed
    rrna_fasta = params.rrna_fasta
    virna_fasta = params.virna_fasta

    // Parse sample sheet
    Channel.fromPath(params.reads)
    | splitCsv(header:true)
    | map { row -> tuple(row.sample_name, row.fastq_read1, row.fastq_read2) }
    | set { read_pairs_ch } // tuple val(name), path(read 1), path(read 2)

    fastqc_ch = read_pairs_ch.flatMap { sample_name, fastq_read1, fastq_read2 ->  
        [[sample_name, fastq_read1], [sample_name, fastq_read2]]  
    }

    // FastQC and MultiQC for raw FASTQs
    FASTQC(fastqc_ch)
    fastqc_out_html = FASTQC.out.html
    fastqc_zip_ch = FASTQC.out.zip
    fastqc_zip_ch.map{ it[1] }.collect()
    | set { multiqc_ch }
    multiqc_report_01 = MULTIQC_01( multiqc_ch )

    // Adapter Trimming
    TRIM_GALORE(read_pairs_ch)
    trimmed_files_ch = TRIM_GALORE.out.trimmed_fastqs
    trim_reports = TRIM_GALORE.out.reports

    // rRNA alignment and viral RNA alignment (sanity check before RiboDetector)
    BOWTIE2_INDEX(rrna_fasta, virna_fasta)
    BOWTIE2_ALIGN(trimmed_files_ch, BOWTIE2_INDEX.out.rrna_index, BOWTIE2_INDEX.out.virna_index)
    rrna_virna_logs = BOWTIE2_ALIGN.out.logs

    // Get viral RNA counts
    VIRAL_COUNTS(BOWTIE2_ALIGN.out.virna_bam) // tuple val(sample_name), path(bam) 
    VIRAL_COUNTS.out.counts.map { name, tsv -> tsv }.collect()
    | set { viral_counts_ch }
    viral_counts_matrix = VIRAL_MATRIX(viral_counts_ch)
    
    // Ribodetector; ML-based ribosomal RNA detection, output filtered fastq's
    // RIBODETECTOR(trimmed_files_ch, params.read_length)
    // filtered_fastqs = RIBODETECTOR.out.filtered_fastqs
    // ribodetector_logs = RIBODETECTOR.out.logs

    // Generate STAR Index and align
    STAR_INDEX(ref_genome, ref_gtf)
    STAR_ALIGN(trimmed_files_ch, STAR_INDEX.out.index) // output bam files tuple val(sample_name), path(bam)
    star_align_out = STAR_ALIGN.out.bam
    STAR_ALIGN.out.log.map{ it[1] }.collect() // collect all STAR output logs into list
    | set { star_logs } 
    STAR_ALIGN.out.gene_counts.map{ it[1] }.collect()
    | set { star_gene_counts }
    
    // Post-alignment MultiQC
    trim_reports
        .mix(star_logs, star_gene_counts, rrna_virna_logs)
        .collect()
        | set { post_align_multiqc_ch }
    post_align_multiqc = MULTIQC_02(post_align_multiqc_ch)

    // Index and Deduplication
    // NOTE: STAR --outSAMtype SoortedByCoordinate means that STAR_ALIGN.out.bam is already sorted, SAMTOOLS sort is redundant
    DEDUPLICATE(STAR_ALIGN.out.bam) // note: doesn't actually remove duplicates, just labels them as duplicates
    dedup_ch = DEDUPLICATE.out.dedup_bam   // (sample, dedup.bam, dedup.bam.bai)
    dedup_logs = DEDUPLICATE.out.log.map { it[1] }

    // Coverage Plots
    // input for bam coverage: tuple val(sample_name), path(bam), path(bai)
    BAMCOVERAGE(dedup_ch)
    bigwig_out = BAMCOVERAGE.out.bigwigs
    
    // Quality control
    dedup_bams_ch = dedup_ch.map { it -> tuple(it[0], it[1]) }
    QUALIMAP(dedup_bams_ch, ref_gtf)
    qualimap_stats = QUALIMAP.out.qualimap_dir
    SAMTOOLS_FLAGSTAT(dedup_bams_ch)
    samtools_stats = SAMTOOLS_FLAGSTAT.out.flagstat
    // RSeQC: gene body coverage and read distribution
    RSEQC_RD(dedup_ch, ref_bed)
    RSEQC_GB(dedup_ch, ref_bed)
    // collect for multiqc
    RSEQC_RD.out.readdist.mix(RSEQC_GB.out.genebody).collect()
    | set { rseqc_logs }

    // Counts Extraction
    all_bams_ch = dedup_ch.map { it -> it[1] }.collect()
    FEATURECOUNTS(all_bams_ch, ref_gtf)
    feature_counts = FEATURECOUNTS.out.counts
    feature_counts_summary = FEATURECOUNTS.out.summary

    SALMON_INDEX(ref_transcriptome)
    SALMON_QUANT(read_pairs_ch, SALMON_INDEX.out.index) // transcriptome-based quantification; note that this will have discrepency with FeatureCounts counts
    salmon_quant_dirs = SALMON_QUANT.out.quant_dir // tuple val(sample_name), path(quant_dir)
    SALMON_QUANT.out.quant_dir
        .map { name, dir -> dir }
        .collect()
        | set { salmon_dirs }

    TXIMPORT(salmon_dirs, file(params.tx2gene))
    salmon_counts = TXIMPORT.out.gene_counts
    salmon_counts_txi = TXIMPORT.out.txi

    DESEQ2(feature_counts, salmon_counts_txi)
    deseq2_outdir = DESEQ2.out.outdir

    // MultiQC
    qualimap_stats
        .mix(samtools_stats, dedup_logs, feature_counts_summary, salmon_dirs, rseqc_logs)
        .collect()
        | set { final_multiqc_ch }
    final_multiqc_report = MULTIQC_03(final_multiqc_ch)

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
    salmon_quant_dirs = salmon_quant_dirs
    salmon_counts = salmon_counts
    salmon_counts_txi = salmon_counts_txi
    deseq2_outdir = deseq2_outdir
    final_multiqc_report = final_multiqc_report
}

output {

    fastqc_out_html {

        path { sample_name, fastqc_html -> "pre-alignment_qc/${sample_name}" }

    }

    fastqc_out_zip {

        path { sample_name, fastqc_zip -> "pre-alignment_qc/${sample_name}" }

    }

    multiqc_report_01 {

        path { multiqc_html -> "pre-alignment_qc" }

    }

    trim_galore_out {

        path { sample_name, read1_trimmed, read2_trimmed -> "trim_galore/${sample_name}" }

    }

    rrna_virna_logs {

        path { logs -> "rRNA_virna_logs"}

    }

    star_align_out {

        path { sample_name, align_files -> "star_output/${sample_name}" }

    }

    post_align_multiqc {

        path { multiqc_html -> "post-alignment_QC" }

    }

    viral_counts_matrix {

        path { viral_counts_matrix -> "viral_counts" }

    }


    bigwig_out {

        path { sample_name, bigwig_fwd, bigwig_rev -> "coverage_plots/${sample_name}" }

    }

    // htseq_counts {

    //     path { sample_name, counts -> "htseq_counts/${sample_name}"}
        
    // }

    feature_counts {

        path { counts -> "feature_counts" }

    }

    // stringtie_counts {

    //     path { counts -> "stringtie_counts" }
    // }

    salmon_quant_dirs {

        path { sample_name, quant_dir -> "salmon_quant/${sample_name}" }

    }

    salmon_counts {

        path { salmon_counts -> "salmon_counts"}

    }

    salmon_counts_txi {

        path { salmon_counts_txi -> "salmon_counts"}

    }

    deseq2_outdir {

        path { deseq2_outdir -> "deseq_outputs" }
        
    }

    final_multiqc_report {

        path { multiqc_html -> "final_QC" }
        
    }

}