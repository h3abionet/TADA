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
