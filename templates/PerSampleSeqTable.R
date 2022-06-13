#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

combineFiles <- list.files(path = '.', pattern = '.${readmode}.RDS\$')
pairIds <- sub('.${readmode}.RDS', '', combineFiles)

combined <- lapply(combineFiles, function (x) readRDS(x))
names(combined) <- pairIds
seqtab <- makeSequenceTable(combined)
seqtab <- seqtab[,nchar(colnames(seqtab)) >= ${params.minLen}]

saveRDS(seqtab, "seqtab.original.${readmode}.RDS")

# this is an optional filtering step to remove *merged* sequences based on 
# min/max length criteria
if (${params.minMergedLen} > 0) {
   seqtab <- seqtab[,nchar(colnames(seqtab)) >= ${params.minMergedLen}]
}

if (${params.maxMergedLen} > 0) {
   seqtab <- seqtab[,nchar(colnames(seqtab)) <= ${params.maxMergedLen}]
}

saveRDS(seqtab, "seqtab.${readmode}.RDS")
saveRDS(combined, "all.${readmode}.RDS")
