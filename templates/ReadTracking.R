#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(dplyr))

getN <- function(x) sum(getUniques(x))

# the gsub here might be a bit brittle...
dadaFs <- as.data.frame(sapply(readRDS("all.dd.R1.RDS"), getN))
rownames(dadaFs) <- gsub('.R1.filtered.fastq.gz', '',rownames(dadaFs))
colnames(dadaFs) <- c("denoised.R1")
dadaFs\$SampleID <- rownames(dadaFs)

# TODO: needs to be optional if R1 only (SE reads)
dadaRs <- as.data.frame(sapply(readRDS("all.dd.R2.RDS"), getN))
rownames(dadaRs) <- gsub('.R2.filtered.fastq.gz', '',rownames(dadaRs))
colnames(dadaRs) <- c("denoised.R2")
dadaRs\$SampleID <- rownames(dadaRs)

# TODO: needs to be optional if no mergers (SE reads)
all.mergers <- readRDS("${mergers}")
mergers <- as.data.frame(sapply(all.mergers, function(x) sum(getUniques(x %>% filter(accept)))))
rownames(mergers) <- gsub('.R1.filtered.fastq.gz', '',rownames(mergers))
colnames(mergers) <- c("merged")
mergers\$SampleID <- rownames(mergers)

seqtab.nochim <- as.data.frame(rowSums(readRDS("${sTable}")))
rownames(seqtab.nochim) <- gsub('.${seqtype}.filtered.fastq.gz', '',rownames(seqtab.nochim))
colnames(seqtab.nochim) <- c("seqtab.${seqtype}.nochim")
seqtab.nochim\$SampleID <- rownames(seqtab.nochim)

trimmed <- read.csv("${trimmedTable}")

track <- Reduce(function(...) merge(..., by = "SampleID",  all.x=TRUE),  list(trimmed, dadaFs, dadaRs, mergers, seqtab.nochim))
# dropped data in later steps gets converted to NA on the join
# these are effectively 0
track[is.na(track)] <- 0

write.table(track, "all.readtracking.${seqtype}.txt", sep = "\\t", row.names = FALSE)
