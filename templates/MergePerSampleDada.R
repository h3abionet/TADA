#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(openssl))

dadaOpt <- ${dadaOpt}

if (!is.na(dadaOpt)) {
  setDadaOpt(dadaOpt)
  cat("dada Options:\\n",dadaOpt,"\\n")
}

dadaopts <- getDadaOpt()

# we want to standardize these for later changes, so let's generate a 
# simple FASTA file of the priors, using 
generate_priors <- function(x, opts = dadaopts, idtype="md5") {
	st <- makeSequenceTable(x)
	# Moving to using the pseudo priors code from Ben here:
    # https://github.com/benjjneb/dada2/blame/278f5f3ec03a846fe157b283cc08f2dd30430ae0/R/dada.R#L400
	priors <- colnames(st)[colSums(st>0) >= opts\$PSEUDO_PREVALENCE | colSums(st) >= opts\$PSEUDO_ABUNDANCE]
	if (length(priors) > 0) {
		ids <- switch(idtype, simple=paste("priorF_", 1:length(priors), sep = ""),
	                            md5=md5(priors))
		seqs.dna <- ShortRead(sread = DNAStringSet(priors), id = BStringSet(ids))
		return(seqs.dna)
	} else {
		return(NA)
	}
}

# this is necessary for QC, but we also want to do this if we want priors from the run
dadaFs <- lapply(list.files(path = '.', pattern = '.dd.R1.RDS'), function (x) readRDS(x))
dadaRs <- lapply(list.files(path = '.', pattern = '.dd.R2.RDS'), function (x) readRDS(x))
names(dadaFs) <- sub('.dd.R1.RDS', '', list.files('.', pattern = '.dd.R1.RDS'))
saveRDS(dadaFs, "all.dd.R1.RDS")

if (as.logical("${params.generate_priors}")) {
	priorsF <- generate_priors(dadaFs, idtype="${params.idType}")
	if (is.na(priorsF)) {
		message("No priors found for R1!")
		file.create("priors.${params.idType}.R1.fna")
	} else {
		writeFasta(priorsF, file = 'priors.${params.idType}.R1.fna')
	}
	if (length(dadaRs) > 0) {
		priorsR <- generate_priors(dadaRs, idtype="${params.idType}")
		if (is.na(priorsR)) {
			message("No priors found for R2!")
			file.create("priors.${params.idType}.R2.fna")
		} else {
			writeFasta(priorsR, file = "priors.${params.idType}.R2.fna")
		}
	}	
}

if (length(dadaRs) > 0) {
	names(dadaRs) <- sub('.dd.R2.RDS', '', list.files('.', pattern = '.dd.R2.RDS'))
	saveRDS(dadaRs, "all.dd.R2.RDS")
	if (as.logical("${params.generate_priors}")) {
		priorsR <- generate_priors(dadaRs, idtype="${params.idType}")
		writeFasta(priorsR, file = 'priors.${params.idType}.R2.fna')
	}	
}