import glob
import re
import os
import datetime
import pandas as pd
import numpy as np
import sys

samples_file = sys.argv[1]

base_dir = os.path.dirname(os.path.abspath(__file__))

# define paths
rrna_results = os.path.join(base_dir, '2_2_rrna_alignment')

# read samples.txt
with open(samples_file, 'r') as f:
    sample_names = [s.strip() for s in f.readlines()]

total_samples = len(sample_names)
print(f'Total samples: {total_samples}')

rrna_files = glob.glob(f'{rrna_results}/*rrna_align.log')

def check_for_missing():

    if len(rrna_files) < 1:
        print(f'No rRNA alignment reports found.')
    else:
        print(f'Found {len(rrna_files)} rRNA alignment reports.')

    reports = []
    total_reports = len(rrna_files)
    missing_reports = 0

    for f in rrna_files:
        f = f.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/2_2_rrna_alignment/', '', f) # this will need to be changed in order to be reproducible
        sample = re.sub(r'_rrna_align\.log$', '', filename) # extracts sample name
        reports.append(sample)

    for s in sample_names:
        if s not in reports:
            print(f"rRNA alignment report for {s} not found.")
            missing_reports += 1

    if total_samples == total_reports:
        print(f'rRNA alignment reports present: {total_reports}/{total_samples}; missing {missing_reports} reports.')
    else:
        print(f'All rRNA alignment reports completed.')


def check_if_empty():
    total_empty = 0

    for file in rrna_files:
        file = file.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/2_2_rrna_alignment/', '', file)
        sample = re.sub(r'_rrna_align\.log$', '', filename) # extracts sample name
        with open(file, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) < 16:
                print(f"rRNA alignment report for {sample} is empty/incomplete.")
                total_empty += 1
    
    if total_empty > 0:
        print(f'{total_empty} files were empty/incomplete.')
    else:
        print(f'No empty/incomplete files found.')

def get_rRNA_percent():

    rrna_by_sample = {}

    for file in rrna_files:
        file = file.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/2_2_rrna_alignment/', '', file)
        sample = re.sub(r'_rrna_align\.log$', '', filename) # extracts sample name
        with open(file, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) < 15:
                print('rRNA alignment report incomplete')
            else:
                rrna_pct = re.sub(r'% overall alignment rate', '', lines[-1]).strip()
                #print(f'{sample} rRNA % overall alignment rate: {rrna_pct}%')
                rrna_pct = float(rrna_pct)
                rrna_by_sample[sample] = rrna_pct

    rrna_pct_df = pd.DataFrame(list(rrna_by_sample.items()), columns=['Sample', 'rRNA Alignment %'])

    rrna_pct_df1 = rrna_pct_df.copy()
    rrna_pct_df1['rRNA Content Status'] = np.where(rrna_pct_df1['rRNA Alignment %'] >= 15, 'High', 'Normal')
    print(f"Samples with high rRNA content: \n{rrna_pct_df1['Sample'][rrna_pct_df1['rRNA Content Status'] == 'High']}")
    problem_samples = [s for s in rrna_pct_df1['Sample'][rrna_pct_df1['rRNA Content Status'] == 'High']]

    # save as csv
    file_path = rrna_results + '/' + f'rRNA_report' + '.csv'
    rrna_pct_df1.to_csv(file_path, index=False)

def main():
    print('Checking for missing rRNA alignment reports...')
    check_for_missing()

    print(f'Checking for empty files...')
    check_if_empty

    print(f'Extracting rRNA alignment stats...')
    get_rRNA_percent()

    print(f'rRNA Content Analysis completed.')

if __name__ == "__main__":
    main()


        



