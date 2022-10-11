#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

tax <- readRDS("${tax}")
map <- readRDS("${map}")

# Note that we use the old ASV ID for output here
write.table(data.frame('ASVID' = row.names(tax), tax),
    file = 'tax_final.${params.idType}.${seqtype}.txt',
    row.names = FALSE,
    col.names=c('#OTU ID', colnames(tax)), sep = "\t")

tax[is.na(tax)] <- "Unclassified"

taxa_combined <- apply(tax, 1, function(x) paste(x, collapse=";"))
taxa_out <- data.frame(names(taxa_combined), taxa_combined)
colnames(taxa_out) <- c("#OTU ID", "taxonomy")

write.table(data.frame('ASVID' = row.names(tax), tax),
    file = 'tax_final.${params.idType}.${seqtype}.full.txt',
    row.names = FALSE,
    col.names=c('#OTU ID', colnames(tax)), sep = "\t")

if (file.exists('bootstrap_final.RDS')) {
    boots <- readRDS("${bt}")
    write.table(data.frame('ASVID' = row.names(boots), boots),
        file = 'tax_final.bootstraps.${params.idType}.${seqtype}.full.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(boots)), sep = "\t")
}

# Write modified data
saveRDS(tax, "tax_final.${params.idType}.${seqtype}.RDS")