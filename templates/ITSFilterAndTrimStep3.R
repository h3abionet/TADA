#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))

out1 <- readRDS("${trimming}")
out2 <- filterAndTrim(fwd = paste0("${pairId}",".R1.cutadapt.fastq.gz"),
                    filt = paste0("${pairId}", ".R1.filtered.fastq.gz"),
                    rev = paste0("${pairId}",".R2.cutadapt.fastq.gz"),
                    filt.rev = paste0("${pairId}", ".R2.filtered.fastq.gz"),
                    maxEE = c(${params.maxEEFor},${params.maxEERev}),
                    truncLen = c(${params.truncFor},${params.truncRev}),
                    truncQ = ${params.truncQ},
                    maxN = ${params.maxN},
                    rm.phix = as.logical(${params.rmPhiX}),
                    maxLen = ${params.maxLen},
                    minLen = ${params.minLen},
                    compress = TRUE,
                    verbose = TRUE,
                    multithread = ${task.cpus})
#Change input read counts to actual raw read counts
out3 <- cbind(out1, out2)
colnames(out3) <- c('input', 'filterN', 'cutadapt', 'filtered')
write.csv(out3, paste0("${pairId}", ".trimmed.txt"))
