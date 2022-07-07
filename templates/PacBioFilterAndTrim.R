#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))

out1 <- removePrimers("${reads}", 
    paste0("${id}",".noprimer.fastq.gz"), 
    primer.fwd="${params.fwdprimer}", 
    primer.rev=dada2:::rc("${params.revprimer}"), 
    orient=TRUE, 
    compress=TRUE,
    verbose=TRUE)

out2 <- filterAndTrim(fwd = paste0("${id}",".noprimer.fastq.gz"),
                    filt = paste0("${id}",".filtered.fastq.gz"),
                    maxEE = ${params.maxEEFor},
                    maxN = ${params.maxN},
                    maxLen = ${params.maxLen},
                    minLen = ${params.minLen},
                    compress = TRUE,
                    verbose = TRUE,
                    multithread = ${task.cpus})
#Change input read counts to actual raw read counts
# out2[1] <- out1[1]
write.csv(out2, paste0("${id}", ".trimmed.txt"))
