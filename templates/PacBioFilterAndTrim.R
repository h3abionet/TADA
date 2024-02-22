#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))

out2 <- filterAndTrim(fwd = "${reads}",
                    filt = "${meta.id}.R1.filtered.fastq.gz",
                    maxEE = ${params.maxEEFor},
                    maxN = ${params.maxN},
                    maxLen = ${params.maxLen},
                    minLen = ${params.minLen},
                    compress = TRUE,
                    verbose = TRUE,
                    multithread = ${task.cpus})

#Change input read counts to actual raw read counts
write.csv(out2, paste0("${meta.id}", ".trimmed.txt"))
