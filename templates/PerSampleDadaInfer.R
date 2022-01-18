#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
  setDadaOpt(dadaOpt)
  cat("dada Options:\\n",dadaOpt,"\\n")
}

errF <- readRDS("errors.R1.RDS")
errR <- readRDS("errors.R2.RDS")
cat("Processing:", "${pairId}", "\\n")

derepF <- derepFastq("${r1}")
ddF <- dada(derepF, err=errF, multithread=${task.cpus}, pool=as.logical("${params.pool}"))

derepR <- derepFastq("${r2}")
ddR <- dada(derepR, err=errR, multithread=${task.cpus}, pool=as.logical("${params.pool}"))

merger <- mergePairs(ddF, derepF, ddR, derepR,
    returnRejects = TRUE,
    minOverlap = ${params.minOverlap},
    maxMismatch = ${params.maxMismatch},
    trimOverhang = as.logical("${params.trimOverhang}"),
    justConcatenate=as.logical("${params.justConcatenate}")
    )

saveRDS(merger, paste("${pairId}.merged.RDS", sep="."))

saveRDS(ddF, "${pairId}.dd.R1.RDS")
# saveRDS(derepFs, "${pairId}.derepFs.RDS")

saveRDS(ddR, "${pairId}.dd.R2.RDS")
# saveRDS(derepRs, "${pairId}.derepRs.RDS")