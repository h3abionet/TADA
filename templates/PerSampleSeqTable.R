#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

mergerFiles <- list.files(path = '.', pattern = '.*.RDS\$')
pairIds <- sub('.merged.RDS', '', mergerFiles)
mergers <- lapply(mergerFiles, function (x) readRDS(x))
names(mergers) <- pairIds
seqtab <- makeSequenceTable(mergers)
seqtab <- seqtab[,nchar(colnames(seqtab)) >= ${params.minLen}]

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
saveRDS(mergers, "all.mergers.RDS")
