#!/bin/bash -l

#$ -P lau-bumc
#$ -N star_index
#$ -l h_rt=12:00:00
#$ -l mem_per_core=16G
#$ -pe omp 8

module load star/2.7.10b

mkdir -p /projectnb/lau-bumc/emily/star_index_aedes

STAR --runMode genomeGenerate \
  --genomeDir /projectnb/lau-bumc/emily/star_index_aedes \
  --genomeFastaFiles /projectnb/lau-bumc/nclau/MosquitoProject/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5.fa \
  --sjdbGTFfile /projectnb/lau-bumc/nclau/MosquitoProject/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.1.gtf \
  --runThreadN 8