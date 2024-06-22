#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(openssl))

# we want to standardize these for later changes, so let's generate a 
# simple FASTA file of the priors, using 
generate_priors <- function(x, idtype) {
	priors <- makeSequenceTable(x) |> colnames()
	ids <- switch(idtype, simple=paste("priorF", 1:length(priors), sep = ""),
                            md5=md5(priors))
	seqs.dna <- ShortRead(sread = DNAStringSet(priors), id = BStringSet(ids))
	return(seqs.dna)
}

# this is necessary for QC, but we also want to do this if we want priors from the run
dadaFs <- lapply(list.files(path = '.', pattern = '.dd.R1.RDS'), function (x) readRDS(x))
names(dadaFs) <- sub('.dd.R1.RDS', '', list.files('.', pattern = '.dd.R1.RDS'))
saveRDS(dadaFs, "all.dd.R1.RDS")

if (as.logical("${params.generate_priors}")) {
	priorsF <- generate_priors(dadaFs, "${params.idType}")
	writeFasta(priorsF, file = 'priors.${params.idType}.R1.fna')
}

dadaRs <- lapply(list.files(path = '.', pattern = '.dd.R2.RDS'), function (x) readRDS(x))

if (length(dadaRs) > 0) {
	names(dadaRs) <- sub('.dd.R2.RDS', '', list.files('.', pattern = '.dd.R2.RDS'))
	saveRDS(dadaRs, "all.dd.R2.RDS")
	if (as.logical("${params.generate_priors}")) {
		priorsR <- generate_priors(dadaRs, "${params.idType}")
		writeFasta(priorsR, file = 'priors.${params.idType}.R2.fna')
	}	
}

