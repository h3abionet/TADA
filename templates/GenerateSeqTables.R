#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
suppressPackageStartupMessages(library(ShortRead))

seqtab <- readRDS("${st}")

if (as.logical('${params.sampleRegex}' != FALSE )) {
    rownames(seqtab) <- gsub('${params.sampleRegex}', "\\\\1", rownames(seqtab), perl = TRUE)
}

# Generate table output
write.table(data.frame('SampleID' = row.names(seqtab), seqtab),
    file = 'seqtab_final.${seqtype}.txt',
    row.names = FALSE,
    col.names=c('#SampleID', colnames(seqtab)), sep = "\\t")

######################################################################
# Convert to simple table + FASTA, from
# https://github.com/LangilleLab/microbiome_helper/blob/master/convert_dada2_out.R#L69
######################################################################

# Generate OTU table output (rows = samples, cols = ASV)
write.table(data.frame('SampleID' = row.names(seqtab), seqtab),
    file = 'seqtab_final.${params.idType}.${seqtype}.txt',
    row.names = FALSE,
    col.names=c('#SampleID', colnames(seqtab)),
    sep = "\\t")

# Generate OTU table for QIIME2 import (rows = ASVs, cols = samples)
write.table(
    data.frame('Taxa' = colnames(seqtab), t(seqtab), check.names = FALSE),
    file = 'seqtab_final.${params.idType}.${seqtype}.qiime2.txt',
    row.names = FALSE,
    quote=FALSE,
    sep = "\\t")

# Write modified data
saveRDS(seqtab, "seqtab_final.${params.idType}.${seqtype}.RDS")
