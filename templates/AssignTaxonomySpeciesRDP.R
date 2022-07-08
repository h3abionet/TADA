#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

seqtab <- readRDS("${st}")

# Assign taxonomy
tax <- NULL
boots <- NULL

if ( ${params.taxBatch} == 0 ) { # no batch, run normally
    cat("Running all samples\n")
    tax <- assignTaxonomy(seqtab, "${ref}",
                    multithread=${task.cpus},
                    tryRC = TRUE,
                    outputBootstraps = TRUE,
                    minBoot = ${params.minBoot},
                    verbose = TRUE)

    boots <- tax\$boot
    if (is.na(as.logical("${sp}"))) {
        tax <- addSpecies(tax\$tax, "${sp}",
             tryRC = TRUE,
             verbose = TRUE)
    }
} else {
    # see https://github.com/benjjneb/dada2/issues/1429 for this
    to_split <- seq(1, ncol(seqtab), by = ${params.taxBatch})
    to_split2 <- c(to_split[2:length(to_split)]-1, ncol(seqtab))

    for(i in 1:length(to_split)){
        cat(paste("Running all samples from",to_split[i], "to", to_split2[i], "\n"))
        seqtab2 <- seqtab[, to_split[i]:to_split2[i]]
        tax2 <- assignTaxonomy(seqtab2, "${ref}",
                multithread=${task.cpus},
                tryRC = TRUE,
                outputBootstraps = TRUE,
                minBoot = ${params.minBoot},
                verbose = TRUE)

        if (is.null(boots)) {
            boots <- tax2\$boot
        } else {
            boots <- rbind(boots, tax2\$boot)
        }

        if (is.na(as.logical("${sp}"))) {
            tax2 <- addSpecies(tax2\$tax, 
                refFasta = "${sp}", 
                tryRC = TRUE,
                verbose = TRUE)
        }
        if (is.null(tax)) {
            tax <- tax2
        } else {
            tax <- rbind(tax, tax2)
        }
    }
}

rownames(tax) <- colnames(seqtab)

# Write original data
saveRDS(tax, "tax_final.RDS")
saveRDS(boots, "bootstrap_final.RDS")
