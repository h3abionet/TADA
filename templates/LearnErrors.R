#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
    setDadaOpt(dadaOpt)
    cat("dada Options:\n",dadaOpt,"\n")
}

# TODO: perform --reads sanity checking (only allow R1 or R2)

# File parsing (these come from the process input channel)
filts <- list.files('.', pattern=paste0("${readmode}",".filtered.fastq.gz"), full.names = TRUE)
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
set.seed(100)

# Learn forward error rates
err <- learnErrors(filts, multithread=${task.cpus})

# This is a rough correction for NovaSeq binning issues
# See https://github.com/h3abionet/TADA/issues/31, we'll likely
# add alternatives here soon

if (as.logical(${params.qualityBinning}) == TRUE ) {
    print("Running binning correction")
    errs <- t(apply(getErrors(err), 1, function(x) { x[x < x[40]] = x[40]; return(x)} ))
    err\$err_out <- errs
}

pdf(paste0("${readmode}",".err.pdf"))
plotErrors(err, nominalQ=TRUE)
dev.off()

saveRDS(err, paste0("errors.","${readmode}",".RDS"))