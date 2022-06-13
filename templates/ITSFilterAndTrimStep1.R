#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))

#Filter out reads with N's
out1 <- filterAndTrim(fwd = "${reads[0]}",
                    filt = paste0("${meta.id}", ".R1.noN.fastq.gz"),
                    rev = if("${reads[1]}" == "null") NULL else "${reads[1]}",
                    filt.rev = if("${reads[1]}" == "null") NULL else paste0("${meta.id}", ".R2.noN.fastq.gz"),
                    maxN = 0,
                    multithread = ${task.cpus})
FWD.RC <- dada2:::rc("${params.fwdprimer}")
REV.RC <- dada2:::rc("${params.revprimer}")

# this may switch to 'env' in the process at some point: 
# https://www.nextflow.io/docs/latest/process.html?highlight=env#output-env
# untested within R though

forP <- file("forward_rc")
writeLines(FWD.RC, forP)
close(forP)

revP <- file("reverse_rc")
writeLines(REV.RC, revP)
close(revP)

saveRDS(out1, "${meta.id}.out.RDS")
