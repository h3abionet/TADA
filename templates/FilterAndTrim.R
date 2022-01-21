#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(optparse))

out <- filterAndTrim(fwd        = "${reads[0]}",
                    filt        = "${pairId}.R1.filtered.fastq.gz",
                    rev         = "${reads[1]}",
                    filt.rev    = "${pairId}.R2.filtered.fastq.gz",
                    trimLeft    = c(${params.trimFor}, ${params.trimRev}),
                    truncLen    = c(${params.truncFor}, ${params.truncRev}),
                    maxEE       = c(${params.maxEEFor}, ${params.maxEERev}),
                    truncQ      = ${params.truncQ},
                    maxN        = ${params.maxN},
                    rm.phix     = as.logical(${params.rmPhiX}),
                    maxLen      = ${params.maxLen},
                    minLen      = ${params.minLen},
                    compress    = TRUE,
                    verbose     = TRUE,
                    multithread = ${task.cpus}
                    )

colnames(out) <- c('input', 'filtered')

write.csv(out, "${pairId}.trimmed.txt")
