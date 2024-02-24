#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(ape))

tree <- read.tree(file = "${tree}")

midtree <- midpoint(tree)

write.tree(midtree, file = "rooted.${seqtype}.newick")
