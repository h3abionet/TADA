#!/usr/bin/env Rscript
library(phangorn)
library(ape)

tree <- read.tree(file = "${tree}")

midtree <- midpoint(tree)

write.tree(midtree, file = "rooted.newick")
