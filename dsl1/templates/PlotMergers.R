#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

mergers <- readRDS("${mergers}")

# filter out any samples w/o reads (nrows == 0)
all_mergers <- mergers |> discard( function(x) nrow(x) == 0 ) |> bind_rows(.id = 'SampleID')

all_mergers\$SampleID <- factor(gsub(".R1.filtered.fastq.gz", "", all_mergers\$SampleID))

# we only keep those that pass here (accept == TRUE)

gg <- all_mergers |> filter(accept) |> ggplot(aes(x = nmatch, y = SampleID)) + 
    stat_bin2d(binwidth=c(2,1), aes(fill = after_stat(ncount))) + 
    scale_fill_viridis_c() + theme_minimal()

if (${params.minMergedLen} > 0) {
	gg <- gg + geom_vline(xintercept=${params.minMergedLen}, linetype="dashed", color = "blue")
}

if (${params.maxMergedLen} > 0) {
	gg <- gg + geom_vline(xintercept=${params.maxMergedLen}, linetype="dashed", color = "red")
}

ggsave('read-overlap-heatmap.pdf', plot=gg, device = 'pdf', width = 7, height = (0.09 * nlevels(all_mergers\$SampleID)), units = 'in')

# save the plot; we may want to make this dynamic (e.g. plotly)
saveRDS(gg, 'read-overlap-heatmap.RDS')
