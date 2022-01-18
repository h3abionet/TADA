#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

#environment(mergePairsRescue) <- asNamespace('dada2')

filtFs <- list.files('.', pattern="R1.filtered.fastq.gz", full.names = TRUE)
filtRs <- list.files('.', pattern="R2.filtered.fastq.gz", full.names = TRUE)

# read in denoised reads for both
ddFs <- readRDS("all.dd.R1.RDS")
ddRs <- readRDS("all.dd.R2.RDS")

mergers <- mergePairs(ddFs, filtFs, ddRs, filtRs,
    returnRejects = TRUE,
    minOverlap = ${params.minOverlap},
    maxMismatch = ${params.maxMismatch},
    trimOverhang = as.logical(${params.trimOverhang}),
    justConcatenate = as.logical(${params.justConcatenate}),
    verbos = TRUE
    )

# TODO: make this a single item list with ID as the name, this is lost
# further on
saveRDS(mergers, "all.mergers.RDS")

# go ahead and make seqtable
seqtab <- makeSequenceTable(mergers)

saveRDS(seqtab, "seqtab.original.RDS")

# this is an optional filtering step to remove *merged* sequences based on 
# min/max length criteria
if (${params.minMergedLen} > 0) {
   seqtab <- seqtab[,nchar(colnames(seqtab)) >= ${params.minMergedLen}]
}

if (${params.maxMergedLen} > 0) {
   seqtab <- seqtab[,nchar(colnames(seqtab)) <= ${params.maxMergedLen}]
}

saveRDS(seqtab, "seqtab.RDS")
