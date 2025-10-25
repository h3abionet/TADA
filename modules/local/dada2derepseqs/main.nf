process DADA2_DEREP_SEQS {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:docker-DADA-1.36"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.RDS"), emit: derep_rds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    derepsF <- derepFastq(${reads[0]}, n=100000, verbose=TRUE)
    derepsR <- derepFastq(${reads[1]}, n=100000, verbose=TRUE)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    """
}
