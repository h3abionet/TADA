#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

seqtabs <- list.files("*.RDS")

seqtab <- mergeSequenceTables(seqtabs, repeats = 'sum', orderBy = 'abundance')

saveRDS(seqtab, "seqtab.merged.RDS")
