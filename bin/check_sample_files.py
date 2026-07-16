import glob
import re
import os
import argparse
import datetime

parser = argparse.ArgumentParser(description='Script for checking on the contents of RNAseq Nextflow output folders.')
parser.add_argument("--process", type=str, default='all', help='Input which process to check, default all') # optional argument
parser.add_argument('--time', type=bool, help='Choose whether to estimate time remaining for STAR process.') # optional argument

args = parser.parse_args()
print(f'Checking {args.process} process...')

# base directory = wherever the script lives
base_dir = os.path.dirname(os.path.abspath(__file__))

# reference other files relative to it
samples_file = os.path.join(base_dir, 'rna_seq_samples.txt')
star_output  = os.path.join(base_dir, '2_star_mapping_output')
fastqc_output = os.path.join(base_dir, '1_fastqc_and_multiqc_reports')
fastq_dir = os.path.join(base_dir, 'fastq_symlinks/')

# read samples.txt
with open('rna_seq_samples.txt', 'r') as f:
    sample_names = [s.strip() for s in f.readlines()]

total_samples = len(sample_names)
print(f'Total samples: {total_samples}')

star_files = glob.glob(f'{star_output}/*Aligned.sortedByCoord.out.bam')
star_complete = []

# CHECK FASTQC RESULT FILES

def check_fastqc_multiqc_progress():
    fastqc_files = glob.glob(f'{fastqc_output}/*fastqc.html') # list of FASTQC HTML reports present in folder

    if len(fastqc_files) < 1:
        print(f'No FASTQC reports found.')
    else:
        print(f'Found {len(fastqc_files)} FASTQC reports.')

    reports = []
    total_reports = 0
    missing_reports = 0

    for f in fastqc_files:
        f = f.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/1_fastqc_and_multiqc_reports/', '', f)
        report_name = re.sub(r'_fastqc\.html$', '', filename) # extracts sample w/ plus R1/R2
        multiqc_report = 'multiqc_report.html'
        reports.append(report_name)

        for s in sample_names:
            report_ct = 0
            for rep in reports:
                if re.sub(r'_R[12]$', '', rep) == s:
                    report_ct += 1
                    total_reports += 1
            if report_ct != 2:
                #print(f'Missing one or more report files for sample: {s}')
                miss = 2 - report_ct
                missing_reports += miss

    if total_samples*2 != total_reports:
        print(f'FASTQC reports present: {total_reports}/{total_samples*2}; missing {missing_reports} reports.')
    else:
        print(f'All FASTQC reports completed.')
    
    if multiqc_report:
        print(f'MultiQC report generated.')
    else:
        print(f'No MultiQC report found.')

# CHECK STAR RESULT FILES

def check_star_progress():
    for f in star_files:
        f = f.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/2_star_mapping_output/', '', f)
        sample_name = re.sub(r'Aligned\.sortedByCoord\.out\.bam$', '', filename) # extracts sample names (e.g. AeAeg_AGO2_F_ovary_ttlRNA_rep1)
        star_complete.append(sample_name)

    print(f'Samples aligned w/ STAR so far: {len(star_complete)}/{total_samples}')

    for sample in sample_names:
        progress_log = os.path.join(star_output, f'{sample}Log.progress.out')
        if sample not in star_complete and os.path.isfile(progress_log):
            mtime = datetime.datetime.fromtimestamp(os.path.getmtime(progress_log))
            print(f'{sample}: RUNNING (last update: {mtime.strftime("%H:%M:%S")})')

# Estimate STAR runtime for remaining samples

def time_left():
    total_time = 0
    progress_logs = glob.glob(f'{star_output}/*Log.progress.out')

    # get most recently modified
    progress_log = max(progress_logs, key=os.path.getmtime)
    sample = os.path.basename(progress_log).replace('Log.progress.out', '')

    r1 = os.path.join(fastq_dir, f'{sample}_R1.fastq')
    if not os.path.isfile(r1):
        print(f'{sample}: R1 file not found')
    else:
        with open(r1) as f:
            line_count = sum(1 for _ in f)        
    
    # check progress log to get average speed for this mapping
    if os.path.isfile(progress_log):
        with open(progress_log, 'r') as f:
            lines = f.readlines()
            speed_tot = 0
            row_ct = 0

            rows = [l for l in lines if l.strip() and l.startswith('May')]
            for val in rows:
                fields = val.split()
                speed = float(fields[3])
                speed_tot += speed
                row_ct += 1

            avg_speed = speed_tot/row_ct
            print(f'Average alignment speed: {avg_speed}')

    reads = line_count // 4
    est_minutes = round((reads / 1e6) / (avg_speed / 60), 1) 
    print(f'{sample}: {reads:,} reads (~{est_minutes} min at {avg_speed}M/hr)')

if args.process == 'fastqc':
    def main():
        check_fastqc_multiqc_progress()
elif args.process == 'star':
    if args.time == True:
        def main():
            check_star_progress()
            time_left()
    else:
        def main():
            check_star_progress()
else:
    if args.time == True:
        def main():
            check_fastqc_multiqc_progress()
            check_star_progress()  
            time_left()
    def main():
        check_fastqc_multiqc_progress()
        check_star_progress()        

if __name__ == "__main__":
    main()
    