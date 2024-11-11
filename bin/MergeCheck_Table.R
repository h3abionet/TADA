#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

mergesum_files <- list.files(".", 
                        pattern = "*.mergesummary.txt", 
                        full.names = TRUE)

mergesums <- lapply(mergesum_files,
                    function(x) { 
                      tmp <- read_tsv(x, 
                               n_max = 2,
                               col_types = "ic",
                               col_names = c("Count", "Metric"))[,1] %>%
                        select(Count) %>% t() %>% as.data.frame()
                      colnames(tmp) <- c("Total", "Merged")
                      rownames(tmp) <- sub("\\S+/(\\S+).mergesummary.txt", "\\1", x)
                      tmp
                      })

names(mergesums) <- gsub("\\S+/(\\S+).mergesummary.txt", "\\1", mergesum_files)

mergesums <- bind_rows(mergesums) %>% mutate(Percent_Merged = round(100 * (Merged/Total), digits = 3 ))

write_tsv(mergesums, "all_mergedstats.tsv", col_names = TRUE)
