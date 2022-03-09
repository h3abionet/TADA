#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

# Note this might be a bit brittle, assumes file ext will be ".RDS"!
# This is primarily to ensure that there isn't an empty ASV sequence,
# which seem to come from merging reads and returning failed merges

seqtabs <- sapply(list.files(pattern = "*.RDS"), function(x) { 
      tmp <- readRDS(x)
      tmp <- tmp[,sapply(colnames(tmp), nchar) > 0]
      return(tmp)
      })

seqtab <- mergeSequenceTables(tables = seqtabs, repeats = 'sum', orderBy = 'abundance')

saveRDS(seqtab, "seqtab.merged.RDS")
