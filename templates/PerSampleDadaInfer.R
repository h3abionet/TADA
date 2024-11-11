#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(Biostrings))

getPriors <- function(x) {
  priors <- readDNAStringSet(x) |> as.vector() |> unname()
  return(priors)
}

set.seed(100)

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
  setDadaOpt(dadaOpt)
  cat("dada Options:\\n",dadaOpt,"\\n")
}

cat("Processing:", "${meta.id}", "\\n")

errF <- readRDS("errors.R1.RDS")
derepF <- derepFastq("${reads[0]}", n=100000)

# TODO: there is probably a better way of doing this 
# when using optparse
paramsF <- list(
  derep=derepF, 
  err=errF,
  multithread=${task.cpus}, 
  pool=as.logical("${params.pool}")
  )
if (as.logical("$run_fpriors")) {
  paramsF\$priors <- getPriors("${fp}")
}

ddF <- do.call(dada, paramsF)
saveRDS(ddF, "${meta.id}.dd.R1.RDS")

if (file.exists ("errors.R2.RDS")) {
  errR <- readRDS("errors.R2.RDS")
  derepR <- derepFastq("${reads[1]}", n=100000)
  paramsR <- list(
    derep=derepR, 
    err=errR, 
    multithread=${task.cpus}, 
    pool=as.logical("${params.pool}")
    )
  if (as.logical("$run_rpriors")) {
    paramsR\$priors <- getPriors("${rp}")
  }

  message("DADA2 params, R2:", paramsR, "\\n")
  ddR <- do.call(dada, paramsR)
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


