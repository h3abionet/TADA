#!/usr/bin/env Rscript
library(dada2)
packageVersion("dada2")

seqtab <- readRDS("${st}")

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "${ref}",
                multithread=${task.cpus},
                tryRC = TRUE,
                outputBootstraps = TRUE,
                minBoot = ${params.minBoot},
                verbose = TRUE)
boots <- tax\$boot

tax <- addSpecies(tax\$tax, "${sp}",
         tryRC = TRUE,
         verbose = TRUE)

rownames(tax) <- colnames(seqtab)

# Write original data
saveRDS(tax, "tax_final.RDS")
saveRDS(boots, "bootstrap_final.RDS")
