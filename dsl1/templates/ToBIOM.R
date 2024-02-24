#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(biomformat))
packageVersion("biomformat")
seqtab <- readRDS("${sTable}")
taxtab <- readRDS("${tTable}")
st.biom <- make_biom(t(seqtab), observation_metadata = taxtab)
write_biom(st.biom, "dada2.${seqtype}.biom")
