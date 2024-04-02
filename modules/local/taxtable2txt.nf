process TAXTABLE2TEXT {

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(tax)
    path(bt)
    path(map)

    output:
    path("tax_final.RDS"), emit: taxtabRDS
    path("tax_final.txt"), emit: taxtab2qiime
    path("*.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(ShortRead))

    tax <- readRDS("${tax}")
    map <- readRDS("${map}")

    # Note that we use the old ASV ID for output here
    write.table(data.frame('ASVID' = row.names(tax), tax),
        file = 'tax_final.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(tax)), sep = "\t")

    tax[is.na(tax)] <- "Unclassified"

    taxa_combined <- apply(tax, 1, function(x) paste(x, collapse=";"))
    taxa_out <- data.frame(names(taxa_combined), taxa_combined)
    colnames(taxa_out) <- c("#OTU ID", "taxonomy")

    write.table(data.frame('ASVID' = row.names(tax), tax),
        file = 'tax_final.full.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(tax)), sep = "\t")

    if (file.exists('bootstrap_final.RDS')) {
        boots <- readRDS("${bt}")
        write.table(data.frame('ASVID' = row.names(boots), boots),
            file = 'tax_final.bootstraps.full.txt',
            row.names = FALSE,
            col.names=c('#OTU ID', colnames(boots)), sep = "\t")
    }

    # Write modified data
    saveRDS(tax, "tax_final.RDS")
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
