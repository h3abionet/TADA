#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))
st.all <- readRDS("${st}")

# Remove chimeras
seqtab <- removeBimeraDenovo(
    st.all, 
    method="consensus", 
    multithread=${task.cpus}, 
    verbose=TRUE ${chimOpts} 
    )

saveRDS(seqtab, "seqtab_final.${seqtype}.RDS")
