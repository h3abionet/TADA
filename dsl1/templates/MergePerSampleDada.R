#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaFs <- lapply(list.files(path = '.', pattern = '.dd.R1.RDS'), function (x) readRDS(x))
names(dadaFs) <- sub('.dd.R1.RDS', '', list.files('.', pattern = '.dd.R1.RDS'))
saveRDS(dadaFs, "all.dd.R1.RDS")

dadaRs <- lapply(list.files(path = '.', pattern = '.dd.R2.RDS'), function (x) readRDS(x))

if (length(dadaRs) > 0) {
	names(dadaRs) <- sub('.dd.R2.RDS', '', list.files('.', pattern = '.dd.R2.RDS'))
	saveRDS(dadaRs, "all.dd.R2.RDS")
}
