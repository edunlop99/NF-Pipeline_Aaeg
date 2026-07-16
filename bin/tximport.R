#!/usr/bin/env Rscript
library(tximport)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
tx2gene_file <- args[1]

tx2gene <- read_tsv(tx2gene_file, col_names = c("TXNAME", "GENEID"))

# quant.sf files are staged into the work dir (current dir)
files <- list.files(".", pattern = "quant.sf$", recursive = TRUE, full.names = TRUE)
names(files) <- sub("_salmon.*", "", basename(dirname(files)))
stopifnot(length(files) > 0, all(file.exists(files)))

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE,
                countsFromAbundance = "no")

write.table(txi$counts, "gene_counts.tsv", sep = "\t", quote = FALSE, col.names = NA)
saveRDS(txi, "txi.rds")