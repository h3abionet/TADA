process PLOTQUALITYPROFILE {
    tag "$meta.id"
    label 'process_low'

    // TODO: pin to a versioned docker instance
    container "ghcr.io/h3abionet/tada:dev"

    input:
    tuple val(meta), path(reads)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.RDS"), emit: rds

    // TODO nf-core: List additional required output channels/values here
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(tidyverse))

    fns <- list.files(pattern="fastq.gz", full.names=TRUE)

    pl <- plotQualityProfile(fns)
    ggsave("${meta.id}.qualities.pdf", plot=pl, device="pdf")

    # we may revisit the quality scores and other info in this plot for other purposes
    saveRDS(pl, "${meta.id}.qualities.RDS")
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.PDF
    """
}
