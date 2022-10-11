#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

seqs <- readRDS("${st}")

# Assign taxonomy
tax <- assignTaxonomy(seqs\$seq, "${ref}",
                multithread=${task.cpus},
                tryRC = TRUE,
                outputBootstraps = TRUE,
                minBoot = ${params.minBoot},
                verbose = TRUE)

boots <- tax\$boot

rownames(tax\$tax) <- seqs[rownames(tax),]\$id
rownames(boots) <- seqs[rownames(boots),]\$id

# Write to disk
saveRDS(tax\$tax, "tax_final.${seqtype}.RDS")
saveRDS(boots, "bootstrap_final.${seqtype}.RDS")
