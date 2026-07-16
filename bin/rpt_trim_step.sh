#!/bin/bash -l

#$ -P lau-bumc
#$ -N trim_repeat
#$ -l h_rt=12:00:00
#$ -l mem_per_core=16G
#$ -pe omp 8

set -o pipefail

set +e
module load miniconda
set -e
conda activate /projectnb/lau-bumc/emily/.conda/envs/rnaseq_pipeline

trim_galore --cores 8 \
            --paired \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGA \
            --output_dir /projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/fastq_trimmed \
            /projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/fastq_symlinks/AeAeg_WT_M_carc_ttlRNA_rep1_R1.fastq \
            /projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/fastq_symlinks/AeAeg_WT_M_carc_ttlRNA_rep1_R2.fastq