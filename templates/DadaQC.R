#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))

fns <- list.files("./fastq", full.names=TRUE)

# this may switch to 'env' in the process at some point: 
# https://www.nextflow.io/docs/latest/process.html?highlight=env#output-env
# untested within R though

pl <- plotQualityProfile(fns, aggregate=TRUE)
ggsave("aggregate-qualities.${readtype}.pdf", plot=pl, device="pdf")

# we may revisit the quality scores and other info in this plot for other purposes
saveRDS(pl, "aggregate-qualities.${readtype}.RDS")
