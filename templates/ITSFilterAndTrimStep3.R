#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))

out1 <- readRDS("${trimming}")
out2 <- filterAndTrim(fwd = paste0("${meta.id}",".R1.cutadapt.fastq.gz"),
                    filt = paste0("${meta.id}", ".R1.filtered.fastq.gz"),
                    rev = if("${reads[1]}" == "null") NULL else paste0("${meta.id}",".R2.cutadapt.fastq.gz"),
                    filt.rev = if("${reads[1]}" == "null") NULL else paste0("${meta.id}", ".R2.filtered.fastq.gz"),
                    maxEE = if("${reads[1]}" == "null") ${params.maxEEFor} else c(${params.maxEEFor}, ${params.maxEERev}), 
                    truncLen = if("${reads[1]}" == "null") ${params.truncFor} else c(${params.truncFor}, ${params.truncRev}),
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
write.csv(out3, paste0("${meta.id}", ".trimmed.txt"))
