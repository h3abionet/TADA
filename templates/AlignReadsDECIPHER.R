#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(DECIPHER))

seqs <- readDNAStringSet("${seqs}")
alignment <- AlignSeqs(seqs,
           anchor=NA,
           processors = ${task.cpus})
writeXStringSet(alignment, "aligned_seqs.fasta")
