#!/usr/bin/env Rscript
library(dada2)
packageVersion("dada2")
library(dplyr)

getN <- function(x) sum(getUniques(x))

# the gsub here might be a bit brittle...
dadaFs <- as.data.frame(sapply(readRDS("all.dd.R1.RDS"), getN))
rownames(dadaFs) <- gsub('.R1.filtered.fastq.gz', '',rownames(dadaFs))
colnames(dadaFs) <- c("denoisedF")
dadaFs\$SampleID <- rownames(dadaFs)

dadaRs <- as.data.frame(sapply(readRDS("all.dd.R2.RDS"), getN))
rownames(dadaRs) <- gsub('.R2.filtered.fastq.gz', '',rownames(dadaRs))
colnames(dadaRs) <- c("denoisedR")
dadaRs\$SampleID <- rownames(dadaRs)

all.mergers <- readRDS("${mergers}")
mergers <- as.data.frame(sapply(all.mergers, function(x) sum(getUniques(x %>% filter(accept)))))
rownames(mergers) <- gsub('.R1.filtered.fastq.gz', '',rownames(mergers))
colnames(mergers) <- c("merged")
mergers\$SampleID <- rownames(mergers)

seqtab.nochim <- as.data.frame(rowSums(readRDS("${sTable}")))
rownames(seqtab.nochim) <- gsub('.R1.filtered.fastq.gz', '',rownames(seqtab.nochim))
colnames(seqtab.nochim) <- c("seqtab.nochim")
seqtab.nochim\$SampleID <- rownames(seqtab.nochim)

trimmed <- read.csv("${trimmedTable}")

track <- Reduce(function(...) merge(..., by = "SampleID",  all.x=TRUE),  list(trimmed, dadaFs, dadaRs, mergers, seqtab.nochim))
# dropped data in later steps gets converted to NA on the join
# these are effectively 0
track[is.na(track)] <- 0

write.table(track, "all.readtracking.txt", sep = "\\t", row.names = FALSE)
