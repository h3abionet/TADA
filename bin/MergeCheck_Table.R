#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

log_files <- list.files(".", 
                        pattern = "*.log", 
                        full.names = TRUE)

process_logs <- function(x) {
    log_lines <- grep("^\\s+", read_lines(x), value = TRUE) 
    log_lines <- gsub("^\\s+|\\s+$", "", log_lines)
    # get rid of parentheses
    log_lines <- gsub("\\s+\\(\\S+\\)$", "", log_lines)
    # at the moment we only keep the paired/merged reads
    merge_data <- str_split_fixed(log_lines, "\\s+", 2) %>% as.data.frame()
    merge_data[,1] <- as.numeric(merge_data[,1])
    read_data <- as.data.frame(t(merge_data[1:2,1]))
    colnames(read_data) <- c("Total", "Merged")
    read_data
}

log_data <- lapply(log_files, process_logs)

names(log_data) <- gsub("\\S+/(\\S+).log", "\\1", log_files)

mergesums <- bind_rows(log_data, .id = "Sample") %>% mutate(Percent_Merged = round(100 * (Merged/Total), digits = 3 ))

write_tsv(mergesums, "all_mergedstats.tsv", col_names = TRUE)