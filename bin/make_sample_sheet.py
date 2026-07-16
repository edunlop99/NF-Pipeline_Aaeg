import glob
import re
import os
import sys
import pandas as pd

# want to make a traditional sample sheet with sample,read1_fastq,read2_fastq' format

fastq_dir = sys.argv[1]

fastq_files = sorted(glob.glob(f"{fastq_dir}/*fastq.gz"))

# Separate R1 and R2
r1_files = [f for f in fastq_files if re.search(r"_R1_001\.(fastq(\.gz)?|fq)$", f)]
r2_files = [f for f in fastq_files if re.search(r"_R2_001\.(fastq(\.gz)?|fq)$", f)]

sample_names = [re.sub(r"_R1_001\.(fastq(\.gz)?|fq)$", "", os.path.basename(f)) for f in r1_files]

# make sample_sheet with sample_name,fastq_r1,fastq_r2 format

# Build samplesheet
samplesheet = pd.DataFrame({
    "sample_name": sample_names,
    "fastq_read1": r1_files,
    "fastq_read2": r2_files
})

samplesheet.to_csv("samplesheet.csv", index=False)