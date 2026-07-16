#!/bin/bash

set -e
set -o pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'  # No color

if [ "$#" -ne 6 ]; then
    echo -e "${RED}Error: Invalid number of arguments.${NC}"
    echo "Usage: $0 <config_directory> <ref_fasta> <fastq_input_path> <read1_suffix> <read2_suffix> <num_threads>"
    exit 1
fi

# Input variables
config_directory=$1
ref_fasta=$2
fastq_input_path=$3
read1_suffix=$4
read2_suffix=$5
num_threads=$6
output_dir=${config_directory}/2_3_vrna_alignment

log_file="${config_directory}/0_nextflow_logs/vrna_alignment.log"
log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "$log_file"
}

if [ ! -d "$output_dir" ]; then
    log "Creating output directory: $output_dir"
    mkdir -p "$output_dir"
fi

# Check if the fastq input path exists
if [ ! -d "$fastq_input_path" ]; then
    echo -e "${RED}Error: Input directory $fastq_input_path does not exist.${NC}"
    exit 1
fi

if [ -z "$num_threads" ]; then
    num_threads=$(nproc)
    log "Number of threads not provided. Using default: $num_threads"
else
    log "Using $num_threads threads."
fi

if [ ! -f "$ref_fasta" ]; then
    echo -e "${RED}Error: Reference FASTA $ref_fasta does not exist.${NC}"
    exit 1
fi

# check for bowtie2 index; if it doesn't exist, create it:
if [ ! -f "$output_dir/vrna_index.1.bt2" ]; then
    log "Bowtie2 vRNA index not found; generating now ..."
    bowtie2-build --threads ${num_threads} ${ref_fasta} ${config_directory}/2_3_vrna_alignment/vrna_index
fi

start_time=$(date +%s)

# Create vRNA counts summary table as .csv
echo -e "sample,virus,length,mapped,unmapped" > ${config_directory}/2_3_vrna_alignment/vrna_counts_summary.csv

# Bowtie2 alignment on all fastq files (not fastq.gz)
# NOTE: output fastq files are named exactly the same as the rRNA filtering output, so it will overwrite them and leave us with double filtered fastq files
log "Starting viral RNA alignment process on files in $fastq_input_path..."
for file in "${fastq_input_path}"/*_filtered${read1_suffix}.fastq; do
    if [[ -f "$file" ]]; then
        sample_name=$(basename "${file%_filtered${read1_suffix}.fastq}")
        R1="$file"
        R2="${file/_filtered${read1_suffix}.fastq/_filtered${read2_suffix}.fastq}"

        expected_report="${config_directory}/2_3_vrna_alignment/${sample_name}_vrna_align.log"
        minimum_size=600
        if [[ -f "$expected_report" ]]; then
            actual_size=$(wc -c <$expected_report)
            if [[ $actual_size -ge $minimum_size ]]; then
                echo "Skipping ${sample_name} — already aligned" | tee -a ${config_directory}/0_nextflow_logs/vrna_alignment.log
                continue
            fi
        fi
        
        log "Processing ${sample_name}..."
        bowtie2 -x ${config_directory}/2_3_vrna_alignment/vrna_index \
                -1 ${R1} \
                -2 ${R2} \
                -p ${num_threads} \
                --un-conc ${config_directory}/fastq_filtered/tmp_${sample_name}_filtered_R%.fastq \
                2> ${config_directory}/2_3_vrna_alignment/${sample_name}_vrna_align.log | \
                samtools view -bS -F 4 | \
                samtools sort -o ${config_directory}/2_3_vrna_alignment/${sample_name}_vrna.sorted.bam
        exit_code=${PIPESTATUS[0]}
        if [ $exit_code -ne 0 ]; then
            rm -f ${config_directory}/fastq_filtered/tmp_${sample_name}_filtered_R*.fastq
            echo "[ERROR] bowtie2 failed for ${sample_name}" | tee -a ${log_file}
            exit 1
        fi

        # overwrite fastq_filtered only if it worked
        mv ${config_directory}/fastq_filtered/tmp_${sample_name}_filtered_R1.fastq \
        ${config_directory}/fastq_filtered/${sample_name}_filtered_R1.fastq
        mv ${config_directory}/fastq_filtered/tmp_${sample_name}_filtered_R2.fastq \
        ${config_directory}/fastq_filtered/${sample_name}_filtered_R2.fastq
        rate=$(grep "overall alignment rate" ${config_directory}/2_3_vrna_alignment/${sample_name}_vrna_align.log | awk '{print $1}')
        echo "Overall vRNA alignment rate for ${sample_name}: ${rate}" | tee -a ${config_directory}/0_nextflow_logs/vrna_alignment.log
        
        # Index and append counts directly to summary
        samtools index ${config_directory}/2_3_vrna_alignment/${sample_name}_vrna.sorted.bam \
        || { echo "[ERROR] samtools index failed for ${sample_name}" | tee -a ${log_file}; exit 1; }
        samtools idxstats ${config_directory}/2_3_vrna_alignment/${sample_name}_vrna.sorted.bam | \
            while IFS=$'\t' read -r virus length mapped unmapped; do
                echo -e "${sample_name},${virus},${length},${mapped},${unmapped}"
            done >> ${config_directory}/2_3_vrna_alignment/vrna_counts_summary.csv

        # Optionally clean up BAM to save storage space
        rm -f ${config_directory}/2_3_vrna_alignment/${sample_name}_vrna.sorted.bam
        rm -f ${config_directory}/2_3_vrna_alignment/${sample_name}_vrna.sorted.bam.bai
    else
        echo -e "${RED}No fastq files found in the directory.${NC}"
        exit 1
    fi
done

log "Successfully completed viral RNA alignment."