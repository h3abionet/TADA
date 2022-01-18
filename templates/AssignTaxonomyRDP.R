#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

seqtab <- readRDS("${st}")

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "${ref}",
                      multithread=${task.cpus},
                      minBoot = ${params.minBoot},
                      tryRC = TRUE,
                      outputBootstraps = TRUE, ${taxLevels}
                      verbose = TRUE 
                      )

# Write to disk
saveRDS(tax\$tax, "tax_final.RDS")
saveRDS(tax\$boot, "bootstrap_final.RDS")
