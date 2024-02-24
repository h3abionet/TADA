#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(DECIPHER))
packageVersion("DECIPHER")

seqtab <- readRDS("${st}")

# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab))

# load database; this should be a RData file
load("${refFile}")

ids <- IdTaxa(dna, trainingSet,
    strand="both",
    processors=${task.cpus},
    verbose=TRUE)
# ranks of interest
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
saveRDS(ids, 'raw_idtaxa.RDS')

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x\$rank)
        taxa <- x\$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab)

boots <- t(sapply(ids, function(x) {
        m <- match(ranks, x\$rank)
        bs <- x\$confidence[m]
        bs
}))
colnames(boots) <- ranks
rownames(boots) <- getSequences(seqtab)

# Write to disk
saveRDS(taxid, "tax_final.${seqtype}.RDS")
saveRDS(boots, "bootstrap_final.${seqtype}.RDS")
