process PER_SAMPLE_SEQTABLE {
    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(mr)
    val(readmode)

    output:
    path("seqtab.${readmode}.RDS"), emit: filtered_seqtable
    path("all.merged.RDS"), optional: true, emit: merged_seqs
    path("seqtab.original.${readmode}.RDS"), emit: seqtabQC
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    combineFiles <- list.files(path = '.', pattern = '.${readmode}.RDS\$')
    pairIds <- sub('.${readmode}.RDS', '', combineFiles)

    combined <- lapply(combineFiles, function (x) readRDS(x))
    names(combined) <- pairIds
    seqtab <- makeSequenceTable(combined)

    saveRDS(seqtab, "seqtab.original.${readmode}.RDS")

    # this is an optional filtering step to remove *merged* sequences based on 
    # min/max length criteria
    if (${params.min_asv_len} > 0) {
       seqtab <- seqtab[,nchar(colnames(seqtab)) >= ${params.min_asv_len}, drop = FALSE]
    }

    if (${params.max_asv_len} > 0) {
       seqtab <- seqtab[,nchar(colnames(seqtab)) <= ${params.min_asv_len}, drop = FALSE]
    }

    saveRDS(seqtab, "seqtab.${readmode}.RDS")
    saveRDS(combined, "all.${readmode}.RDS")
    """
}
