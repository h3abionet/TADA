#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

# note this has to be a seqtable with the actual sequences in them 
# (e.g. IDs are not transformed to md5 or otherwise)
seqtab <- readRDS("${seqtab}")

seqlens <- data.frame(seqs = colnames(seqtab), lengths = nchar(colnames(seqtab)))

gg <- ggplot(seqlens, aes(x = lengths)) + 
	geom_density() + 
	ggtitle("Sequence Length Distribution") + 
	xlab("Length (nt)")

ggsave('asv-length-distribution.pdf', device = 'pdf', height = 3, width = 5, units = 'in')

# save the plot; we may want to make this dynamic (e.g. plotly)
saveRDS(gg, 'asv-length-distribution.RDS')