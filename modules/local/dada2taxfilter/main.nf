process DADA2_TAXFILTER {
    tag "tax_filter:${params.tax_filter_rank}"
    label 'process_single'

    container container "ghcr.io/h3abionet/tada:dev"

    input:
    path readmap
    path seqtab
    path taxtab
    path asvs

    output:
    path "readmap.tax_filtered.tsv", emit: readmap_tax_filtered
    path "seqtab_final.tax_filtered.tsv", emit: seqtab_tax_filtered
    path "taxtab_final.tax_filtered.tsv", emit: taxtab_tax_filtered
    path "asvs.tax_filtered.tsv", emit: asvs_tax_filtered
    path "*.tax_filtered.RDS", emit: rdata
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def rank = params.tax_filter_rank
    """
    #!/usr/bin/env Rscript

    # TODO: quick hack script, could be cleaned up
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(ShortRead))
    suppressPackageStartupMessages(library(tidyverse))

    seqtab <- readRDS("${seqtab}")
    taxtab <- readRDS("${taxtab}")
    readmap <- readRDS("${readmap}")
    taxtab_filtered <- taxtab %>% as_tibble(rownames = "TaxID") %>% filter(Phylum != "Unclassified")
    ids <- taxtab_filtered\$TaxID
    seqtab_filtered <- seqtab[,ids]
    # asvs_filtered <- asvs[ids]
    readmap_filtered <- readmap[readmap$id %in% ids,]
    asvs_filtered <- DNAStringSet(readmap_filtered$seq)
    names(asvs_filtered) <- readmap_filtered$id

    # Generate table output
    write.table(data.frame('SampleID' = row.names(seqtab_filtered), seqtab_filtered),
        file = "seqtab_final.tax_filtered.tsv",
        row.names = FALSE,
        col.names=c('#SampleID', colnames(seqtab_filtered)), sep = "\\t")
    
    writeXStringSet(asvs_filtered, file="asvs.tax_filtered.fna")

    write.table(taxtab_filtered,
        file = "taxtab_final.tax_filtered.tsv",
        row.names = FALSE,
        col.names=colnames(taxtab_filtered), 
        sep = "\\t")

    # Write modified data
    saveRDS(seqtab_filtered, "seqtab_final.tax_filtered.RDS")
    saveRDS(taxtab_filtered, "taxtab_final.tax_filtered.RDS")
    saveRDS(readmap_filtered, "readmap_final.tax_filtered.RDS")
    """

    // stub:
    // def args = task.ext.args ?: ''
    
    // // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    // //               Have a look at the following examples:
    // //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    // //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // """
    // touch ${prefix}.bam

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     dada2taxfilter: \$(samtools --version |& sed '1!d ; s/samtools //')
    // END_VERSIONS
    // """
}
