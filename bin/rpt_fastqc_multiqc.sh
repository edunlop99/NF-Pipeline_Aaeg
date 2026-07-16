#!/bin/bash -l

#$ -P lau-bumc
#$ -N qc_repeat
#$ -l h_rt=48:00:00
#$ -l mem_per_core=16G
#$ -pe omp 8

set -e
set -o pipefail

module load fastqc
module load multiqc

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Input variables
config_directory="/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline"
samples=("AeAeg_DCR2_M_head_ttlRNA_rep2" "AeAeg_DCR2_M_testes_ttlRNA_rep2" "AeAeg_WT_M_testes_ttlRNA_rep2" "AeAeg_AGO2_M_testes_ttlRNA_rep2" "AeAeg_R2D2_F_head_ttlRNA_rep2" "AeAeg_AGO2_F_ovary_ttlRNA_rep2")
fastq_input_path="${config_directory}/2_2_rrna_alignment"
num_threads=$NSLOTS
output_dir="${config_directory}/rpt_qc_reports"
map_output_fpath="${config_directory}/2_star_mapping_output_rpt"
read1_suffix="_R1"
read2_suffix="_R2"
log_file="${config_directory}/0_nextflow_logs/rpt_qc_reports.log"

log() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1" | tee -a "${log_file}"
}

log_error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $1${NC}" | tee -a "${log_file}"
}

trap 'log_error "Failed on line $LINENO. Exiting..."' ERR

mkdir -p "${output_dir}"

# Check if the fastq input path exists
if [ ! -d "$fastq_input_path" ]; then
    log_error "Input directory $fastq_input_path does not exist."
    exit 1
fi

if [ -z "$num_threads" ]; then
    num_threads=$(nproc)
    log "Number of threads not provided. Using default: $num_threads"
else
    log "Using $num_threads threads."
fi

start_time=$(date +%s)

# Run FastQC on all fastq files
log "Running FastQC on files in $fastq_input_path..."
for file in "${fastq_input_path}"/*.fastq; do
    if [[ -f "$file" ]]; then
        filename=$(basename "$file")
        expected_report="${output_dir}/${filename%.fastq}_fastqc.html"

        if [ -f "$expected_report" ]; then
            log "Skipping $filename — FastQC report already exists"
            continue
        fi

        log "Processing $filename..."
        fastqc -t $num_threads "$file" --outdir="$output_dir" > /dev/null 2>> ${log_file} \
            || { log_error "FastQC failed for $filename"; exit 1; }
    else
        log_error "No fastq files found in $fastq_input_path"
        exit 1
    fi
done

log "Successfully completed FastQC reports."

# Run MultiQC on FastQC reports
log "Running MultiQC in $output_dir..."
multiqc "${output_dir}" -o "${output_dir}" > /dev/null 2>> ${log_file} \
    || { log_error "MultiQC failed"; exit 1; }

log "Successfully completed MultiQC reports."

end_time=$(date +%s)
execution_time=$((end_time - start_time))

log "${GREEN}FastQC and MultiQC reports generated successfully.${NC}"
log "${GREEN}Reports are stored at: $output_dir${NC}"
log "Total execution time: $(($execution_time / 60)) minutes and $(($execution_time % 60)) seconds."

# Post-alignment mapping metrics
log "Starting mapping metrics generation at $(date)"

if [ ! -f "${config_directory}/generate_map_metrics.py" ]; then
    log_error "Mapping metrics script not found: ${config_directory}/generate_map_metrics.py"
    exit 1
fi

if [ ! -d "${map_output_fpath}" ]; then
    log_error "STAR mapping output directory not found: ${map_output_fpath}"
    exit 1
fi

log "Running mapping metrics:"
log "  Config directory: ${config_directory}"
log "  Mapped files directory: ${map_output_fpath}"
log "python ${config_directory}/generate_map_metrics.py ${config_directory} ${map_output_fpath}"

python "${config_directory}/generate_map_metrics.py" "${config_directory}" "${map_output_fpath}" > /dev/null 2>> ${log_file} \
    || { log_error "generate_map_metrics.py failed"; exit 1; }

log "Running MultiQC on mapping output..."
multiqc "${map_output_fpath}" -o "${map_output_fpath}" > /dev/null 2>> ${log_file} \
    || { log_error "MultiQC on mapping output failed"; exit 1; }

mkdir -p "${output_dir}/map_metrics"
cp -r "${map_output_fpath}/multiqc_report.html" "${map_output_fpath}/multiqc_data" "${output_dir}/map_metrics/" \
    || { log_error "Failed to copy MultiQC reports"; exit 1; }
rm -rf "${map_output_fpath}/multiqc_report.html" "${map_output_fpath}/multiqc_data" 2>> ${log_file}

log "Mapping metrics completed. Reports stored at: ${output_dir}/map_metrics"
log "All QC steps completed successfully."