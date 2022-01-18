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

derepF <- derepFastq("${filtFor}")

ddF <- dada(fileF, err=errF, multithread=${task.cpus}, pool=as.logical("${params.pool}"))

derepR <- derepFastq("${filtRev}")
ddR <- dada(fileR, err=errR, multithread=${task.cpus}, pool=as.logical("${params.pool}"))

merger <- mergePairs(ddF, derepF, ddR, derepR,
    returnRejects = TRUE,
    minOverlap = ${params.minOverlap},
    maxMismatch = ${params.maxMismatch},
    trimOverhang = as.logical("${params.trimOverhang}"),
    justConcatenate=as.logical("${params.justConcatenate}")
    )

saveRDS(merger, paste("${pairId}.merged.RDS", sep="."))

# saveRDS(ddFs, "${pairId}.dd.R1.RDS")
# saveRDS(derepFs, "${pairId}.derepFs.RDS")

# saveRDS(ddRs, "${pairId}.dd.R2.RDS")
# saveRDS(derepRs, "${pairId}.derepRs.RDS")