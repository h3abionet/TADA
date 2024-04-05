#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option_list = list(
    make_option(c("-i","--id"),
        type="character", 
        help="ID of sample"),
    make_option(c("-f", "--file_pattern"),
        type="character",
        default = "fastq.gz",
        help="File pattern [default \"%default\"]"),
    make_option(c("--session"),
        action = "store_true", 
        help="Generate nf-core compliant YAML file w/ version information",
        ),
    make_option(c("--yaml"),
        default = "versions.yml", 
        help="YAML file name (see --session)")
)

opt <- parse_args(OptionParser(option_list=option_list))

fns <- list.files(pattern=opt$file_pattern, full.names=TRUE)

pl <- plotQualityProfile(fns)
ggsave(paste0(opt$id,".qualities.pdf"), plot=pl, device="pdf")

# we may revisit the quality scores and other info in this plot for other purposes
saveRDS(pl, paste0(opt$id,".qualities.RDS"))