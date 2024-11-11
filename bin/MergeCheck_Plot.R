#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

# A very simple plot script for generating a 'heatmap' checking overlaps. 
option_list = list(
    # make_option(c("--forward_clip"), type="character", default=NULL, help="Forward (5') trim"),
    # make_option(c("--reverse_cli["), type="character", default=NULL, help="Reverse (5') trim"),
    make_option(c("--minMergedLen"), type="numeric", default=0, help="cpus"),
)

opt <- parse_args(OptionParser(option_list=option_list))minMergedLen
len_files <- list.files(".", 
                        pattern = "*.lengthstats.txt", 
                        full.names = TRUE)

lens_tmp <- lapply(len_files,
                   read_tsv, 
                   col_names = c("Length", "Count"), col_types = "ii")

names(lens_tmp) <- gsub("\\S+/(\\S+).lengthstats.txt", "\\1", len_files)

# bind all the data, then group by Sample, add in relative abundance per Sample and binning info, then group by Sample + Bin and summarize counts and RelAb per bin
# It's a bit of a hack but it generally works; however it's not perfect
lens_all <- bind_rows(lens_tmp, .id = "Sample") %>% 
  group_by(Sample) %>% 
  mutate(RelAb=Count/sum(Count),
         Bin = cut(Length, 
                   seq(min(Length), 
                       max(Length),5), 
                   include.lowest = TRUE)) %>% 
  group_by(Sample, Bin) %>% 
  mutate(ReadCountPerBin=sum(Count),
         RelAbPerBin=sum(RelAb)) 

# This will become settable, but essentially anything 50nt or less is not kept
cutoff <- 50

gg <- lens_all |> ggplot(aes(x=Length, y=Sample, fill=ReadCountPerBin)) + 
  geom_tile(stat = "identity") +
  scale_fill_viridis_c(option="plasma", direction = -1) +
  annotate("rect",
           xmin = 0,
           xmax = cutoff,
           ymin = 0.5,
           ymax = Inf,
           alpha=0.3, fill="blue") +
  theme_minimal()

# if(fprimer > 0) {
#   gg <- gg+ geom_vline(xintercept=fprimer, color = "red", alpha = 0.5)
# }

# if(rprimer > 0) {
#   gg <- gg+ geom_vline(xintercept=rprimer, color = "black", alpha = 0.5)
# }

# if(maxsizeprimers > 0) {
#   gg <- gg+ geom_vline(xintercept=maxsizeprimers, color = "green", alpha = 0.5)
# }

ggsave(gg, "MergedCheck_heatmap.pdf")
saveRDS(gg, "MergedCheck_heatmap.RDS")
