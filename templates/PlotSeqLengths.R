#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

# change to input channel name
seqtab <- readRDS('seqtab.original.merged.RDS')

seqlens <- data.frame(seqs = colnames(seqtab), lengths = nchar(colnames(seqtab)))

ggplot(seqlens, aes(x = lengths)) + 
	geom_density() + 
	ggtitle("Sequence Length Distribution") + 
	xlab("Length (nt)")

ggsave('asv-length-distribution.png', device = 'png', height = 3, width = 5, units = 'in')
