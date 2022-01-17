#!/usr/bin/env Rscript
library(dada2); packageVersion("dada2")

fns <- list.files("./fastq", full.names=TRUE)

# this may switch to 'env' in the process at some point: 
# https://www.nextflow.io/docs/latest/process.html?highlight=env#output-env
# untested within R though

pdf("qualities.pdf", onefile = TRUE)
for (i in seq(1, length(fns), by = 4)) {
    pl <- plotQualityProfile(fns[i:(i+3)])
    print(pl)
}

dev.off()