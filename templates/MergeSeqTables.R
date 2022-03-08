#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

# note this might be a bit brittle, assumes file ext will be ".RDS"
seqtabs <- sapply(list.files(pattern = "./*.RDS"), function(x) { 
      tmp <- readRDS(x)
      tmp <- tmp[,sapply(colnames(tmp), nchar) > 0]
      return(tmp)
      })

seqtab <- mergeSequenceTables(tables = seqtabs, repeats = 'sum', orderBy = 'abundance')

saveRDS(seqtab, "seqtab.merged.RDS")
