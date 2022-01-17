#!/usr/bin/env Rscript
library(dada2)
library(ShortRead)

tax <- readRDS("${tax}")
map <- readRDS("${map}")

# Note that we use the old ASV ID for output here
write.table(data.frame('ASVID' = row.names(tax), tax),
    file = 'tax_final.txt',
    row.names = FALSE,
    col.names=c('#OTU ID', colnames(tax)), sep = "\t")

# Tax table
if(!identical(rownames(tax), as.character(map\$seq))){
    stop("sequences in taxa and sequence table are not ordered the same.")
}

tax[is.na(tax)] <- "Unclassified"
rownames(tax) <- map\$id
taxa_combined <- apply(tax, 1, function(x) paste(x, collapse=";"))
taxa_out <- data.frame(names(taxa_combined), taxa_combined)
colnames(taxa_out) <- c("#OTU ID", "taxonomy")

write.table(data.frame('ASVID' = row.names(tax), tax),
    file = 'tax_final.simple.full.txt',
    row.names = FALSE,
    col.names=c('#OTU ID', colnames(tax)), sep = "\t")

write.table(taxa_out,
    file = 'tax_final.simple.txt',
    row.names = FALSE,
    sep = "\t")

if (file.exists('bootstrap_final.RDS')) {
    boots <- readRDS("${bt}")
    if(!identical(rownames(boots), as.character(map\$seq))){
        stop("sequences in bootstrap and sequence table are not ordered the same.")
    }
    rownames(boots) <- map\$id
    write.table(data.frame('ASVID' = row.names(boots), boots),
        file = 'tax_final.bootstraps.simple.full.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(boots)), sep = "\t")
}

# Write modified data
saveRDS(tax, "tax_final.simple.RDS")