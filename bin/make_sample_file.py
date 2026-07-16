import glob
import re
import os

file_path_1 = '/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/rna_seq_samples.txt'
file_path_2 = '/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/rna_seq_samples_2.txt'

# Ensure the file exists before checking its size
if os.path.isfile(file_path_1) and os.path.getsize(file_path_1) != 0:
    print(f"{file_path_1} exists and is not empty.")
else:
    print(f"{file_path_1}does not exist or is empty.")
    print("Creating sample file...")

    # all fastq files
    file_paths_rep1 = glob.glob('/projectnb/lau-bumc/nclau/MosquitoProject/AeAeg_KMEL_RNAiMutants_2025/totalRNA_libs_rep1/OrigFastqFiles/*.fastq')
    file_paths_rep2 = glob.glob('/projectnb/lau-bumc/nclau/MosquitoProject/AeAeg_KMEL_RNAiMutants_2025/totalRNA_libs_rep2/OrigFastqFiles/*.fastq')

    samples = []
    for file in file_paths_rep1:
        file = file.strip()
        filename = re.sub("/projectnb/lau-bumc/nclau/MosquitoProject/AeAeg_KMEL_RNAiMutants_2025/totalRNA_libs_rep1/OrigFastqFiles/", "", file)
        sample = re.sub(r'_R[12]\.fastq(\.gz)?$', '', filename) # removes either _R1 or _R2 for each sample
        samples.append(sample)

    for file in file_paths_rep2:
        file = file.strip()
        filename = re.sub(r'/projectnb/lau-bumc/nclau/MosquitoProject/AeAeg_KMEL_RNAiMutants_2025/totalRNA_libs_rep[12]/OrigFastqFiles/', '', file)
        sample = re.sub(r'_R[12]\.fastq(\.gz)?$', '', filename) # removes either _R1 or _R2 for each sample
        samples.append(sample)

    samples = set(samples) # remove duplicates

    with open('/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/rna_seq_samples.txt', 'w') as f:
        for sample in samples:
            f.write(f'{sample}\n')
    
    if os.path.isfile(file_path_1) and os.path.getsize(file_path_1) != 0:
        print(f"{file_path_1} created successfully.")
    else:
        print(f"ERROR: {file_path_1} not created succesfully.")

# Ensure the file exists before checking its size
if os.path.isfile(file_path_2) and os.path.getsize(file_path_2) != 0:
    print(f"{file_path_2} exists and is not empty.")
else:
    print(f"{file_path_2} does not exist or is empty.")
    print("Creating sample file...")

    # all fastq files
    file_paths = glob.glob('/projectnb/lau-bumc/nclau/MosquitoProject/AeAeg_mRNAseq_ZZ20220112/*.fastq')

    samples = []
    for file in file_paths:
        file = file.strip()
        filename = re.sub("/projectnb/lau-bumc/nclau/MosquitoProject/AeAeg_mRNAseq_ZZ20220112/", "", file)
        sample = re.sub(r'_R[12]\.fastq(\.gz)?$', '', filename) # removes either _R1 or _R2 for each sample
        samples.append(sample)

    samples = set(samples) # remove duplicates

    with open('/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/rna_seq_samples_2.txt', 'w') as f:
        for sample in samples:
            f.write(f'{sample}\n')

    if os.path.isfile(file_path_2) and os.path.getsize(file_path_2) != 0:
        print(f"{file_path_2} created successfully.")
    else:
        print(f"ERROR: {file_path_2} not created succesfully.")
        