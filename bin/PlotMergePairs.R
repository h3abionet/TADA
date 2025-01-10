#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# A very simple plot script for generating a 'heatmap' checking overlaps. 
option_list = list(
    make_option(c("--forward_primer"), type="character", default=NULL, help="Forward primer sequence"),
    make_option(c("--reverse_primer"), type="character", default=NULL, help="Reverse primer sequence"),
    make_option(c("--minOverlap"), type="numeric", default=0, help="cpus"),
)

opt <- parse_args(OptionParser(option_list=option_list))

# the files here come from a simple awk one-liner:
# awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, "\t", lengths[l]}}' > SAMPLE.lengthstats.txt
len_files <- list.files(".", pattern = "*.lengthstats.txt", full.names = TRUE)

lens_tmp <- lapply(len_files, read_tsv, col_names = c("Length", "Count"), col_types = cols())
names(lens_tmp) <- gsub("\\S+/(\\S+).lengthstats.txt", "\\1", len_files)
lens_all <- bind_rows(lens_tmp, .id = "Sample")
head(lens_all)

lens_all |> ggplot(aes(x = Length, y = Sample)) + 
  stat_bin2d(binwidth=c(2,1), aes(fill = after_stat(ncount))) + 
  scale_fill_viridis_c(option = "turbo") +
  annotate("rect", 
           xmin = 1, 
           xmax = fprimer + rprimer + cutoff, 
           ymin = 0,
           ymax = Inf,
           alpha=0.4, fill="red") + 
  theme_minimal()


# eval(parse(text=paste0("setDadaOpt(", opt$dadaOptStr, ")")))
# 
# environment(mergePairsRescue) <- asNamespace('dada2')
# 
# filtFs <- list.files('.', pattern="R1.filtered.fastq.gz", full.names = TRUE)
# filtRs <- list.files('.', pattern="R2.filtered.fastq.gz", full.names = TRUE)
# 
# errF <- readRDS(opt$errFor)
# errR <- readRDS(opt$errRev)
# cat("Processing all samples\n")
# 
# #Variable selection from CLI input flag --pool
# pool <- opt$pool
# if(pool == "T" || pool == "TRUE"){
#   pool <- as.logical(pool)
# }
# 
# derepFs <- derepFastq(filtFs)
# 
# ddFs <- dada(derepFs, err=errF, multithread=opt$cpus, pool=pool)
# 
# derepRs <- derepFastq(filtRs)
# 
# ddRs <- dada(derepRs, err=errR, multithread=opt$cpus, pool=pool)
# 
# mergers <- mergePairs(ddFs, derepFs, ddRs, derepRs,
#     returnRejects = TRUE,
#     minOverlap = opt$minOverlap,
#     maxMismatch = opt$maxMismatch,
#     trimOverhang = as.logical(opt$trimOverhang),
#     justConcatenate = as.logical(opt$justConcatenate)
#     )
# 
# # TODO: make this a single item list with ID as the name, this is lost
# # further on
# saveRDS(mergers, "all.mergers.RDS")
# 
# saveRDS(ddFs, "all.ddF.RDS")
# saveRDS(derepFs, "all.derepFs.RDS")
# 
# saveRDS(ddRs, "all.ddR.RDS")
# saveRDS(derepRs, "all.derepRs.RDS")
