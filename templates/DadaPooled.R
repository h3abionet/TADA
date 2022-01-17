#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
  setDadaOpt(dadaOpt)
  cat("dada Options:\\n",dadaOpt,"\\n")
}

filts <- list.files('.', pattern="${readmode}.filtered.fastq.gz", full.names = TRUE)

err <- readRDS("${err}")
cat("Processing all samples\\n")

#Variable selection from CLI input flag --pool
pool <- "${params.pool}"

# 'pool' is a weird flag, either 'pseudo' (string), or T/F (bool)
if(pool != "pseudo"){
  pool <- as.logical(pool)
}

cat(paste0("Denoising ${readmode} reads: pool:", pool, "\\n"))
dds <- dada(filts, err=err, multithread=${task.cpus}, pool=pool)

saveRDS(dds, "all.dd.${readmode}.RDS")
