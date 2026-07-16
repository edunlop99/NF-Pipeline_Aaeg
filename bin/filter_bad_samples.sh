#!/bin/bash -l

#$ -P lau-bumc
#$ -N star_index_rpt
#$ -l h_rt=48:00:00
#$ -l mem_per_core=16G
#$ -pe omp 8

set -e
set -o pipefail

module load star/2.7.10b
module load bowtie2/2.5.1

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

config_directory="/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline"
samples=("AeAeg_DCR2_M_head_ttlRNA_rep2" "AeAeg_DCR2_M_testes_ttlRNA_rep2" "AeAeg_WT_M_testes_ttlRNA_rep2" "AeAeg_AGO2_M_testes_ttlRNA_rep2" "AeAeg_R2D2_F_head_ttlRNA_rep2" "AeAeg_AGO2_F_ovary_ttlRNA_rep2")
fastq_files="${config_directory}/fastq_symlinks"
read1_suffix="_R1"
read2_suffix="_R2"
num_threads=$NSLOTS
star_index_path="/projectnb/lau-bumc/emily/star_index_aedes"
clip5_num=11
clip3_num=5
map_output_fpath="${config_directory}/2_star_mapping_output_rpt"
map_log_fpath="${map_output_fpath}/logs"
log_file="${config_directory}/0_nextflow_logs/rpt_mapping.log"

log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "${log_file}"
}

log_error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $1${NC}" | tee -a "${log_file}"
}

trap 'log_error "Failed on line $LINENO. Exiting..."' ERR

mkdir -p "${map_output_fpath}" "${map_log_fpath}"

log "Starting rRNA filtering and remapping for ${#samples[@]} samples"

for sample_name in "${samples[@]}"; do
    log "Processing rRNA alignment for ${sample_name}..."
    bowtie2 -x ${config_directory}/2_2_rrna_alignment/rrna_index \
        -1 ${fastq_files}/${sample_name}${read1_suffix}.fastq \
        -2 ${fastq_files}/${sample_name}${read2_suffix}.fastq \
        -p ${num_threads} \
        --un-conc ${config_directory}/2_2_rrna_alignment/${sample_name}_filtered_R%.fastq \
        -S /dev/null \
        2> ${config_directory}/2_2_rrna_alignment/${sample_name}_rrna_align.log \
        || { log_error "bowtie2 failed for ${sample_name}"; exit 1; }

    rate=$(grep "overall alignment rate" ${config_directory}/2_2_rrna_alignment/${sample_name}_rrna_align.log | awk '{print $1}')
    log "rRNA alignment rate for ${sample_name}: ${rate}"

    if [ ! -f "${config_directory}/2_2_rrna_alignment/${sample_name}_filtered_R1.fastq" ]; then
        log_error "Filtered fastq not created for ${sample_name}"
        exit 1
    fi

    log "Re-Aligning ${sample_name} with STAR..."
    star_cmd="STAR --runThreadN ${num_threads} \
        --runMode alignReads \
        --genomeDir ${star_index_path} \
        --readFilesIn ${config_directory}/2_2_rrna_alignment/${sample_name}_filtered_R1.fastq \
                    ${config_directory}/2_2_rrna_alignment/${sample_name}_filtered_R2.fastq \
        --clip5pNbases ${clip5_num} ${clip5_num} \
        --clip3pNbases ${clip3_num} ${clip3_num} \
        --outFileNamePrefix ${map_output_fpath}/${sample_name} \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts"

    log "STAR command: ${star_cmd}"

    eval ${star_cmd} >> ${map_log_fpath}/${sample_name}_rpt_mapping_step.log 2>&1 \
        || { log_error "STAR failed for ${sample_name} - check ${map_log_fpath}/${sample_name}_rpt_mapping_step.log"; exit 1; }

    log "Remapping completed for ${sample_name}"

done

log "All samples processed successfully"