process READ_TRACKING {
    label 'process_low'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(tracking)

    output:
    path("readtracking.txt"), emit: read_tracking

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(dplyr))

    getN <- function(x) sum(getUniques(x))

    # this is the main DADA2 trimmed table 
    # (consolidated from filterAndTrim)
    track <- read.csv("all.trimmed.csv")

    # the gsub here might be a bit brittle...
    dadaFs <- as.data.frame(sapply(readRDS("all.dd.R1.RDS"), getN))
    rownames(dadaFs) <- gsub('.R1.filtered.fastq.gz', '',rownames(dadaFs))
    colnames(dadaFs) <- c("denoised.R1")
    dadaFs\$SampleID <- rownames(dadaFs)
    track <- merge(track, dadaFs, by = "SampleID",  all.x=TRUE)

    # TODO: needs to be optional if R1 only (SE reads)
    if ( file.exists("all.dd.R2.RDS") ) {
        dadaRs <- as.data.frame(sapply(readRDS("all.dd.R2.RDS"), getN))
        rownames(dadaRs) <- gsub('.R2.filtered.fastq.gz', '',rownames(dadaRs))
        colnames(dadaRs) <- c("denoised.R2")
        dadaRs\$SampleID <- rownames(dadaRs)
        track <- merge(track, dadaRs, by = "SampleID",  all.x=TRUE)
    }

    # TODO: needs to be optional if no merged data (SE reads)
    if (file.exists("all.merged.RDS")) {
        all.mergers <- readRDS("all.merged.RDS")
        mergers <- as.data.frame(sapply(all.mergers, function(x) sum(getUniques(x %>% filter(accept)))))
        rownames(mergers) <- gsub('.R1.filtered.fastq.gz', '',rownames(mergers))
        colnames(mergers) <- c("merged")
        mergers\$SampleID <- rownames(mergers)
        track <- merge(track, mergers, by = "SampleID",  all.x=TRUE)
    }

    seqtab.nochim <- as.data.frame(rowSums(readRDS("seqtab.${params.id_type}.RDS")))
    rownames(seqtab.nochim) <- gsub('.R1.filtered.fastq.gz', '',rownames(seqtab.nochim))
    colnames(seqtab.nochim) <- c("seqtab.nochim")
    seqtab.nochim\$SampleID <- rownames(seqtab.nochim)
    track <- merge(track, seqtab.nochim, by = "SampleID",  all.x=TRUE)

    # TODO: the below are repetitive with the above, should make this a function
    if (file.exists("seqtab.${params.id_type}.search_filtered.RDS")) {
        seqtab.sf <- as.data.frame(rowSums(readRDS("seqtab.${params.id_type}.search_filtered.RDS")))
        rownames(seqtab.sf) <- gsub('.R1.filtered.fastq.gz', '',rownames(seqtab.sf))
        colnames(seqtab.sf) <- c("seqtab.searchfilter")
        seqtab.sf\$SampleID <- rownames(seqtab.sf)
        track <- merge(track, seqtab.sf, by = "SampleID",  all.x=TRUE)
    }

    if (file.exists("seqtab.tax_filtered.RDS")) {
        seqtab.tf <- as.data.frame(rowSums(readRDS("seqtab.tax_filtered.RDS")))
        rownames(seqtab.tf) <- gsub('.R1.filtered.fastq.gz', '',rownames(seqtab.tf))
        colnames(seqtab.tf) <- c("seqtab.taxfilter")
        seqtab.tf\$SampleID <- rownames(seqtab.tf)
        track <- merge(track, seqtab.tf, by = "SampleID",  all.x=TRUE)
    }

    # dropped data in later steps gets converted to NA on the join
    # these are effectively 0
    track[is.na(track)] <- 0

    write.table(track, "readtracking.txt", sep = "\\t", row.names = FALSE)
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
