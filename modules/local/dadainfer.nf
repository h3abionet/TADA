// TODO: rename file dada2pooledinfer?
process DADA2_POOLED_INFER {
    tag "${readmode}: ${params.pool}"
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:docker-DADA-1.36"

    input:
    tuple val(readmode), file(err), file(reads)

    output:
    path("all.dd.*.RDS"), emit: inferred

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bandsize = params.platform == 'pacbio' ? ', BAND_SIZE=32' : ''
    def dadaOpt = params.dada_opts ? "${params.dada_opts}" : ""
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    dadaOpt <- ${dadaOpt}

    if (!is.na(dadaOpt)) {
      setDadaOpt(dadaOpt)
      cat("dada Options:\\n",dadaOpt,"\\n")
    }
    set.seed(100)

    filts <- list.files('.', pattern="${readmode}.filtered.fastq.gz", full.names = TRUE)

    err <- readRDS("${err}")
    cat("Processing all samples\\n")

    #Variable selection from CLI input flag --pool
    pool <- "${params.pool}"

    # 'pool' is a weird flag, either 'pseudo' (string), or T/F (bool)
    if(pool != "pseudo"){
      pool <- as.logical(pool)
    }

    dereps <- derepFastq(filts, n=100000, verbose=TRUE)

    cat(paste0("Denoising ${readmode} reads: pool:", pool, "\\n"))
    dds <- dada(dereps, 
      err=err, 
      multithread=${task.cpus}, 
      pool=pool ${bandsize})

    saveRDS(dds, "all.dd.${readmode}.RDS")
    """

    stub:
    def args = task.ext.args ?: ''
    """
    # add some real stuff here
    touch all.dd.${readmode}.RDS
    """
}
