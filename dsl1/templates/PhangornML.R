#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(phangorn))

phang.align <- read.phyDat("${aln}",
                            format = "fasta",
                            type = "DNA")

dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
write.tree(fit\$tree, file = "tree.${seqtype}.newick")

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR, "phangorn.tree.${seqtype}.RDS")
write.tree(fitGTR\$tree, file = "tree.GTR.${seqtype}.newick")