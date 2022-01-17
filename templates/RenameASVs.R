#!/usr/bin/env Rscript
library(dada2)
library(ShortRead)
library(digest)

# read RDS w/ data
st <- readRDS("${st}")
st.raw <- readRDS("${rawst}")

# get sequences
seqs <- colnames(st)
seqs.raw <- colnames(st.raw)

# get IDs based on idType
ids_study <- switch("${params.idType}", simple=paste("ASV", 1:ncol(st), sep = ""),
                            md5=sapply(colnames(st), digest, algo="md5"))
ids_study.raw <- switch("${params.idType}", simple=paste("ASV", 1:ncol(st.raw), sep = ""),
                            md5=sapply(colnames(st.raw), digest, algo="md5"))

# sub IDs
colnames(st) <- ids_study
colnames(st.raw) <- ids_study.raw

# generate FASTA
seqs.dna <- ShortRead(sread = DNAStringSet(seqs), id = BStringSet(ids_study))
# Write out fasta file.
writeFasta(seqs.dna, file = 'asvs.${params.idType}.nochim.fna')

seqs.dna.raw <- ShortRead(sread = DNAStringSet(seqs.raw), id = BStringSet(ids_study.raw))
writeFasta(seqs.dna.raw, file = 'asvs.${params.idType}.raw.fna')

# Write modified data (note we only keep the no-chimera reads for the next stage)
saveRDS(st, "seqtab_final.simple.RDS")
saveRDS(data.frame(id = ids_study, seq = seqs), "readmap.RDS")