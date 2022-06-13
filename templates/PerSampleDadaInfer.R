#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
  setDadaOpt(dadaOpt)
  cat("dada Options:\\n",dadaOpt,"\\n")
}

cat("Processing:", "${meta.id}", "\\n")

errF <- readRDS("errors.R1.RDS")
derepF <- derepFastq("${reads[0]}")
ddF <- dada(derepF, err=errF, multithread=${task.cpus}, pool=as.logical("${params.pool}"))
saveRDS(ddF, "${meta.id}.dd.R1.RDS")

if (file.exists ("errors.R2.RDS")) {
  errR <- readRDS("errors.R2.RDS")
  derepR <- derepFastq("${reads[1]}")
  ddR <- dada(derepR, err=errR, multithread=${task.cpus}, pool=as.logical("${params.pool}"))
  saveRDS(ddR, "${meta.id}.dd.R2.RDS")

  merger <- mergePairs(ddF, derepF, ddR, derepR,
      returnRejects = TRUE,
      minOverlap = ${params.minOverlap},
      maxMismatch = ${params.maxMismatch},
      trimOverhang = as.logical("${params.trimOverhang}"),
      justConcatenate=as.logical("${params.justConcatenate}")
      )

  saveRDS(merger, paste("${meta.id}.merged.RDS", sep="."))
} else {
  # yes this is a little silly (it's the same as the dd.R1.RDS above).
  # But it does make the logical flow through this channel easier
  saveRDS(ddF, paste("${meta.id}.R1.RDS", sep="."))
}


