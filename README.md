# NF-Pipeline_Aaeg
Modular Nextflow pipeline for total RNAseq analysis; developed for research purposes.

Emily Dunlop, Summer 2026 Bioinformatics intern
Lau Laboratory, BUMC 
https://class.bu.edu/~nclau/LauLabBUMC/

HOW TO USE:
1) Create input "samplesheet.csv"
  - header: sample_name,read1_fastq,read2_fastq
  - should contain sample name, path to read1 fastq file (can be .gz or .fastq), and path to read2 fastq file.
2) Obtain the following reference files for the organism:
  - Genome FASTA
  - Transcriptome FASTA
  - GTF
  - rRNA reference FASTA (if desired)
  - viral RNA reference FASTA (if desired)
3) Edit nextflow.config file:
  - enter name for "outputDir" at the top; results files will be published here.
  - params block to point to samplesheet.csv and all necessary reference files
  - profile block to add path to conda environtment, change to custom settings. [CURRENTLY CONFIGURED FOR SGE CLUSTER AND CONDA ENVS ONLY]
4) Run with following command:
  - nextflow run main.nf -c nextflow.config -profile cluster,conda

