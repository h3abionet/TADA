#!/usr/bin/env Rscript

# Main purpose of this script is to merge all trimming data into one table

trimmedFiles <- list.files(path = '.', pattern = '*.trimmed.txt')
sample.names <- sub('.trimmed.txt', '', trimmedFiles)
trimmed <- do.call("rbind", lapply(trimmedFiles, function (x) as.data.frame(read.csv(x))))
colnames(trimmed)[1] <- "Sequence"
trimmed\$SampleID <- sample.names
write.csv(trimmed, "all.trimmed.csv", row.names = FALSE)