#!/bin/bash

# extract transcript_id -> gene_id from the GTF
# edit: removes "unassigned_transcript" lines which Salmon cannot read 
# edit: removes isoform version since the ignoreTxVersion flag was not working
grep -v '^#' /projectnb/lau-bumc/emily/NF-Pipeline/refs/AeAeg_GCF_002204515.2.fixed.gtf \
  | awk -F'\t' '$3=="transcript"' \
  | grep -o 'gene_id "[^"]*"; transcript_id "[^"]*"' \
  | sed 's/gene_id "//; s/"; transcript_id "/\t/; s/"//' \
  | awk -F'\t' '$2 !~ /^unassigned_transcript/ {sub(/\.[0-9]+$/, "", $2); print $2"\t"$1}' \
  > tx2gene.tsv

# check overlap
cut -f1 tx2gene.tsv | sed 's/\.[0-9]*$//' | sort -u > /tmp/tx.txt
cut -f1 /projectnb/lau-bumc/emily/NF-Pipeline/work/b7/ab82bca53e5205044c5ae25567b809/AeAeg_Aag2_siGFP_2_LRNA_salmon/quant.sf | tail -n +2 | sed 's/\.[0-9]*$//' | sort -u > /tmp/sal.txt
comm -12 /tmp/sal.txt /tmp/tx.txt | wc -l
wc -l /tmp/sal.txt /tmp/tx.txt