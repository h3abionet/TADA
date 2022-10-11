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

tax <- addSpecies(tax\$tax, "${sp}",
         tryRC = TRUE,
         verbose = TRUE)

# make sure these are the same order
# they should be, but we don't assume this
rownames(tax) <- seqs[rownames(tax),]\$id
rownames(boots) <- seqs[rownames(boots),]\$id

# Write original data
saveRDS(tax, "tax_final.${seqtype}.RDS")
saveRDS(boots, "bootstrap_final.${seqtype}.RDS")
