process DECIPHER {
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(seqs)

    output:
    path("aligned_seqs.fna"), optional: true, emit: alignment

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(DECIPHER))

    seqs <- readDNAStringSet("${seqs}")
    alignment <- AlignSeqs(seqs,
               anchor=NA,
               processors = ${task.cpus})
    writeXStringSet(alignment, "aligned_seqs.fna")
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
