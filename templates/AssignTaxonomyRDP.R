#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

seqs <- readRDS("${st}")
seqtab <- seqs\$seq

# Assign taxonomy
tax <- NULL
boots <- NULL

if ( ${params.taxBatch} == 0 | length(seqtab) < ${params.taxBatch} ) { # no batch, run normally
    cat("Running all samples\\n")
    tax <- assignTaxonomy(seqtab, "${ref}",
                    multithread=${task.cpus},
                    tryRC = TRUE,
                    outputBootstraps = TRUE,
                    minBoot = ${params.minBoot},
                    verbose = TRUE)

    boots <- tax\$boot
} else {
    # see https://github.com/benjjneb/dada2/issues/1429 for this
    to_split <- seq(1, length(seqtab), by = ${params.taxBatch})
    to_split2 <- c(to_split[2:length(to_split)]-1, length(seqtab))

    for(i in 1:length(to_split)){
        cat(paste("Running all samples from",to_split[i], "to", to_split2[i], "\\n"))
        seqtab2 <- seqtab[to_split[i]:to_split2[i]]
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

        if (is.null(tax)) {
            tax <- tax2
        } else {
            tax <- rbind(tax, tax2)
        }
    }
}

# make sure these are the same order
# they should be, but we don't assume this
rownames(tax) <- seqs[rownames(tax),]\$id
rownames(boots) <- seqs[rownames(boots),]\$id

# Write original data
saveRDS(tax, "tax_final.${seqtype}.RDS")
saveRDS(boots, "bootstrap_final.${seqtype}.RDS")
