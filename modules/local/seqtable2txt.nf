process DADA2_SEQTABLE2TEXT {

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(seqtab)

    output:
    path("seqtab_final.qiime2.txt"), emit: seqtab2qiime
    path("*.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(ShortRead))

    seqtab <- readRDS("${seqtab}")

    # Generate table output
    write.table(data.frame('SampleID' = row.names(seqtab), seqtab),
        file = 'seqtab_final.txt',
        row.names = FALSE,
        col.names=c('#SampleID', colnames(seqtab)), sep = "\\t")

    # Generate OTU table for QIIME2 import (rows = ASVs, cols = samples)
    write.table(
        data.frame('Taxa' = colnames(seqtab), t(seqtab), check.names = FALSE),
        file = 'seqtab_final.qiime2.txt',
        row.names = FALSE,
        quote=FALSE,
        sep = "\\t")

    # Write modified data
    saveRDS(seqtab, "seqtab_final.RDS")

    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
