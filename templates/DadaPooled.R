#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
  setDadaOpt(dadaOpt)
  cat("dada Options:\\n",dadaOpt,"\\n")
}

set.seed(100)

filts <- list.files('.', pattern="${readmode}.filtered.fastq.gz", full.names = TRUE)

err <- readRDS("${err}")
cat("Processing all samples\\n")

#Variable selection from CLI input flag --pool
pool <- "${params.pool}"

# 'pool' is a weird flag, either 'pseudo' (string), or T/F (bool)
if(pool != "pseudo"){
  pool <- as.logical(pool)
}

dereps <- derepFastq(filts, n=100000, verbose=TRUE)

cat(paste0("Denoising ${readmode} reads: pool:", pool, "\\n"))
dds <- dada(dereps, err=err, multithread=${task.cpus}, pool=pool)

if (length(filts) == 1) {
  nm <- basename(filts[1])
  dds <- list(nm = dds)
  names(dds)[1] <- nm
}

saveRDS(dds, "all.dd.${readmode}.RDS")
