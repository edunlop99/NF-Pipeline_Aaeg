#!/usr/bin/env Rscript

library(DESeq2)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

fc_deseq2_results <- args[1]
slmn_deseq2_results <- args[2]

