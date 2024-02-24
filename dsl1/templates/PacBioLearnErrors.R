#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
    setDadaOpt(dadaOpt)
    cat("dada Options:\n",dadaOpt,"\n")
}
 
# File parsing
filts <- list.files('.', pattern="filtered.fastq.gz", full.names = TRUE)
sample.namesF <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
set.seed(100)

dereps <- derepFastq(filts, n=${derepreads}, verbose=TRUE)

# Learn forward error rates
errs <- learnErrors(dereps, 
    nbases = 1e8, 
    errorEstimationFunction=PacBioErrfun, 
    BAND_SIZE=32, 
    multithread=${task.cpus},
    verbose=TRUE)

pdf(paste0("${readmode}",".err.pdf"))
plotErrors(errs, nominalQ=TRUE)
dev.off()

saveRDS(errs, paste0("errors.","${readmode}",".RDS"))
