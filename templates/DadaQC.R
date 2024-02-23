#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))

fns <- list.files(pattern="fastq", full.names=TRUE)

# this may switch to 'env' in the process at some point: 
# https://www.nextflow.io/docs/latest/process.html?highlight=env#output-env
# untested within R though

# pl <- plotQualityProfile(fns, aggregate=TRUE)
# ggsave("${meta.id}.aggregate-qualities.pdf", plot=pl, device="pdf")

# # we may revisit the quality scores and other info in this plot for other purposes
# saveRDS(pl, "${meta.id}.aggregate-qualities.RDS")

pl <- plotQualityProfile(fns)
ggsave("${meta.id}.qualities.pdf", plot=pl, device="pdf")

# we may revisit the quality scores and other info in this plot for other purposes
saveRDS(pl, "${meta.id}.qualities.RDS")
