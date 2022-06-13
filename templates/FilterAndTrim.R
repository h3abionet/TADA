#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

out <- filterAndTrim(fwd        = "${reads[0]}",
                    filt        = "${meta.id}.R1.filtered.fastq.gz",
                    rev         = if("${reads[1]}" == "null") NULL else "${reads[1]}",
                    filt.rev    = if("${reads[1]}" == "null") NULL else "${meta.id}.R2.filtered.fastq.gz",
                    trimLeft    = if("${reads[1]}" == "null") ${params.trimFor} else  c(${params.trimFor}, ${params.trimRev}),
                    truncLen    = if("${reads[1]}" == "null") ${params.truncFor} else c(${params.truncFor}, ${params.truncRev}),
                    maxEE       = if("${reads[1]}" == "null") ${params.maxEEFor} else c(${params.maxEEFor}, ${params.maxEERev}), 
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

write.csv(out, "${meta.id}.trimmed.txt")
