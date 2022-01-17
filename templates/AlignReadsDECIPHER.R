#!/usr/bin/env Rscript
library(dada2)
library(DECIPHER)

seqs <- readDNAStringSet("${seqs}")
alignment <- AlignSeqs(seqs,
           anchor=NA,
           processors = ${task.cpus})
writeXStringSet(alignment, "aligned_seqs.fasta")
