#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
    setDadaOpt(dadaOpt)
    cat("dada Options:\n",dadaOpt,"\n")
}
setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})
 
# File parsing
filts <- list.files('.', pattern="filtered.fastq.gz", full.names = TRUE)
sample.namesF <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
set.seed(100)

# Learn forward error rates
errs <- learnErrors(filts, 
    errorEstimationFunction=PacBioErrfun, 
    BAND_SIZE=32, 
    multithread=${task.cpus})

pdf("PacBio.err.pdf")
plotErrors(errs, nominalQ=TRUE)
dev.off()

saveRDS(errs, "errors.RDS")