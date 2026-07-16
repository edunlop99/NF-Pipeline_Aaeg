#!/bin/bash -l

#$ -P lau-bumc
#$ -N vrna_repeat
#$ -l h_rt=24:00:00
#$ -l mem_per_core=16G
#$ -pe omp 16

set -o pipefail

set +e
module load miniconda
set -e
conda activate /projectnb/lau-bumc/emily/.conda/envs/rnaseq_pipeline

config_directory="/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline"
fastq_files="${config_directory}/fastq_trimmed" # gonna have to use these since the ones in fastq_filtered are already filtered from vRNA oof
vrna_fasta="${config_directory}/refs/combined_viral.fa"
read1_suffix="_R1"
read2_suffix="_R2"
is_paired_end=true
num_threads=$NSLOTS

mkdir -p ${config_directory}/2_3_vrna_alignment_rpt
mkdir -p ${config_directory}/fastq_filtered

log_file="${config_directory}/0_nextflow_logs/vrna_alignment_rpt.log"
log() {
    echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] $1" | tee -a "$log_file"
}

echo "[INFO] Starting viral RNA alignment at $(date)" | tee -a $log_file

if [ ! -f "${vrna_fasta}" ]; then
    echo "[ERROR] Reference FASTA not found: ${vrna_fasta}" | tee -a $log_file
    exit 1
fi

if [ ! -d "${fastq_files}" ]; then
    echo "[ERROR] Input directory not found: ${fastq_files}" | tee -a $log_file
    exit 1
fi

# use index from original directory
if [ ! -f "${config_directory}/2_3_vrna_alignment/vrna_index.1.bt2" ]; then
    log "Bowtie2 vRNA index not found; generating now..."
    log "bowtie2-build --threads ${num_threads} ${vrna_fasta} ${config_directory}/2_3_vrna_alignment/vrna_index"
    bowtie2-build --threads ${num_threads} ${vrna_fasta} ${config_directory}/2_3_vrna_alignment/vrna_index
fi

# will output summary to new output directory
if [ ! -f "${config_directory}/2_3_vrna_alignment_rpt/vrna_counts_summary.csv" ]; then
    echo "sample,virus,length,mapped,unmapped" > ${config_directory}/2_3_vrna_alignment_rpt/vrna_counts_summary.csv
fi

log "Checking RNA alignment data on files in ${fastq_files}..."

# input files have format: "samplename_{R1/R2}_val{1/2}.fq"
for file in ${fastq_files}/*${read1_suffix}_val_1.fq; do
    if [[ -f "$file" ]]; then
        sample_name=$(basename "${file%${read1_suffix}_val_1.fq}")
        R1="$file"
        R2="${file/${read1_suffix}_val_1.fq/${read2_suffix}_val_2.fq}"

        # Skip if already processed
        expected_report="${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna_align.log"
        minimum_size=600
        if [[ -f "$expected_report" ]]; then
            actual_size=$(wc -c < "$expected_report")
            if [[ $actual_size -ge $minimum_size ]]; then
                log "Skipping ${sample_name} — already aligned"
                continue
            fi
        fi

        log "Processing ${sample_name}..."
        bowtie2 -x ${config_directory}/2_3_vrna_alignment/vrna_index \
                -1 $R1 \
                -2 $R2 \
                -p ${num_threads} \
                2> ${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna_align.log | \
                samtools view -bS -F 4 | \
                samtools sort -o ${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna.sorted.bam
        exit_code=${PIPESTATUS[0]}
        if [ $exit_code -ne 0 ]; then
            log "[ERROR] bowtie2 failed for ${sample_name}"
            exit 1
        fi

        rate=$(grep "overall alignment rate" ${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna_align.log | awk '{print $1}')
        log "Overall vRNA alignment rate for ${sample_name}: ${rate}"

        set +e # temp disable error exit
        log "Indexing BAM for ${sample_name}..."
        samtools index ${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna.sorted.bam 
        if [ $? -ne 0 ]; then
            log "[ERROR] samtools index failed for ${sample_name}"
            exit 1
        fi
        log "Running idxstats for ${sample_name}..."
        samtools idxstats ${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna.sorted.bam | \
            grep -v "^\*" | \
            while IFS=$'\t' read -r virus length mapped unmapped; do
                echo "${sample_name},${virus},${length},${mapped},${unmapped}"
            done >> ${config_directory}/2_3_vrna_alignment_rpt/vrna_counts_summary.csv
        if [ ${PIPESTATUS[0]} -ne 0 ]; then
            log "[ERROR] samtools idxstats failed for ${sample_name}"
            exit 1
        fi
        set -e
        log "Clearing BAM files for ${sample_name}..."
        rm -f ${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna.sorted.bam
        rm -f ${config_directory}/2_3_vrna_alignment_rpt/${sample_name}_vrna.sorted.bam.bai
        
    else
        log "[ERROR] No fastq files matching pattern found in ${fastq_files}"
        exit 1
    fi

done

log "Viral RNA alignment completed."