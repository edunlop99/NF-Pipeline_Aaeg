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
vrna_results = os.path.join(base_dir, '2_3_vrna_alignment')
counts_summary = os.path.join(base_dir, '2_3_vrna_alignment/vrna_counts_summary.csv')

# read samples.txt
with open(samples_file, 'r') as f:
    sample_names = [s.strip() for s in f.readlines()]

total_samples = len(sample_names)
print(f'Total samples: {total_samples}')

vrna_results = glob.glob(f'{vrna_results}/*_vrna_align.log')

def check_for_missing():

    if len(vrna_results) < 1:
        print(f'No viral RNA alignment reports found.')
    else:
        print(f'Found {len(vrna_results)} viral RNA alignment reports.')

    reports = []
    total_reports = len(vrna_results)
    missing_reports = 0

    for f in vrna_results:
        f = f.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/2_3_vrna_alignment/', '', f)
        sample = re.sub(r'_vrna_align\.log$', '', filename) # extracts sample name
        reports.append(sample)

    for s in sample_names:
        if s not in reports:
            print(f"viral RNA alignment report for {s} not found.")
            missing_reports += 1

    if total_samples == total_reports:
        print(f'viral RNA alignment reports present: {total_reports}/{total_samples}; missing {missing_reports} reports.')
    else:
        print(f'All viral RNA alignment reports completed.')


def check_if_empty():
    total_empty = 0

    for file in vrna_results:
        file = file.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/2_3_vrna_alignment/', '', file)
        sample = re.sub(r'_vrna_align\.log$', '', filename) # extracts sample name
        with open(file, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) < 15:
                print(f"viral RNA alignment report for {sample} is empty/incomplete.")
                total_empty += 1
    
    if total_empty > 0:
        print(f'{total_empty} files were empty/incomplete.')
    else:
        print(f'No empty/incomplete files found.')

def get_vRNA_percent():

    vrna_by_sample = {}

    for file in vrna_results:
        file = file.strip()
        filename = re.sub(r'^/projectnb/lau-bumc/emily/Bulk-RNA-Seq-Nextflow-Pipeline/2_3_vrna_alignment/', '', file)
        sample = re.sub(r'_vrna_align\.log$', '', filename) # extracts sample name
        with open(file, 'r') as f:
            lines = [s.strip() for s in f.readlines()]
            if len(lines) < 15:
                print('')
            else:
                vrna_pct = re.sub(r'% overall alignment rate', '', lines[-1]).strip()
                #print(f'{sample} viral RNA % overall alignment rate: {vrna_pct}%')
                vrna_pct = float(vrna_pct)
                vrna_by_sample[sample] = vrna_pct

    vrna_pct_df = pd.DataFrame(list(vrna_by_sample.items()), columns=['Sample', 'Overall viral RNA Alignment %'])

    # save as csv
    file_path = vrna_results + '/' + f'vRNA_report' + '.csv'
    vrna_pct_df1.to_csv(file_path, index=False)

def main():
    print('Checking for missing viral RNA alignment reports...')
    check_for_missing()

    print(f'Checking for empty files...')
    check_if_empty

    print(f'Extracting viral RNA alignment stats...')
    get_vRNA_percent()

    print(f'viral RNA Content Analysis completed.')

if __name__ == "__main__":
    main()


        



