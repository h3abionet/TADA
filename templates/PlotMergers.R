#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

# this needs to be added to the Docker image
library(ggridges)

mergers <- readRDS("all.merged.RDS"))

# filter out any samples w/o reads (nrows == 0)
all_mergers <- mergers |> discard( function(x) nrow(x) == 0 ) |> bind_rows(.id = 'SampleID')

all_mergers$SampleID <- factor(gsub(".R1.filtered.fastq.gz", "", all_mergers$SampleID))

# we only keep those that pass here (accept == TRUE)
all_mergers |> filter(accept) |> ggplot(aes(x = nmatch, y = SampleID)) + 
	geom_density_ridges(alpha=0.6, bandwidth=3) + theme(axis.text.y=element_text(size=6)) + ggtitle('Paired read overlap (bases)')

ggsave('overlap-size.png', device = 'png', width = 7, height = (0.09 * nlevels(all_mergers\$SampleID)), units = 'in')
