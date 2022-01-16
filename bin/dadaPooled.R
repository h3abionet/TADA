#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("--err"), type="character", default=NULL, help="Error RDS file"),
    make_option(c("--read"), type="character", default=NULL, help="Read (must be R1 or R2)"),
    make_option(c("--cpus"), type="numeric", , default=1, help="cpus"),
    make_option(c("--pool"), type="character", help="Pooling"),
    make_option(c("--dadaOptStr"), type="character", default='', help="values for setDadaOpt passed as one string")
)

opt <- parse_args(OptionParser(option_list=option_list))
ddOpts <- paste0("setDadaOpt(", opt$dadaOptStr, ")")

#cat("Options:\n", opt, "\n")
cat("dada Options:\n", ddOpts, "\n")

#eval(parse(text=ddOpts))

filts <- list.files('.', pattern=paste0(opt$read,".filtered.fastq.gz"), full.names = TRUE)

err <- readRDS(opt$err)
cat("Processing all samples\n")

#Variable selection from CLI input flag --pool
pool <- opt$pool
if(pool == "T" || pool == "TRUE"){
  pool <- as.logical(pool)
}

cat(paste0("Denoising ",opt$read," reads: pool:", pool, "\n"))
dds <- dada(filts, err=err, multithread=opt$cpus, pool=pool)

saveRDS(dds, paste0("all.dd",ifelse(opt$reads == 'R1', 'F', 'R'),".RDS"))
