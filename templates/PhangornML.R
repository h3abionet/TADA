#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(phangorn))

phang.align <- read.phyDat("aligned_seqs.fasta",
                            format = "fasta",
                            type = "DNA")

dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
write.tree(fit\$tree, file = "tree.newick")

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR, "phangorn.tree.RDS")
write.tree(fitGTR\$tree, file = "tree.GTR.newick")