#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

fns <- list.files("./fastq", full.names=TRUE)

# this may switch to 'env' in the process at some point: 
# https://www.nextflow.io/docs/latest/process.html?highlight=env#output-env
# untested within R though

pdf(paste0("aggregate-qualities.pdf"), onefile = TRUE)
pl <- plotQualityProfile(fns, aggregate=TRUE)
ggsave("aggregate-qualities.pdf", plot=pl, device=".pdf")

# we may revisit the quality scores and other info in this plot for other purposes
saveRDS(pl, "aggregate-qualities.RDS")
