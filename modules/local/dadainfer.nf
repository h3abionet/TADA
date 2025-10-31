process DADA2_POOLED_INFER {
    tag "${readmode}: ${params.pool}"
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:docker-DADA-1.36"

    input:
    tuple val(readmode), path(err), path(dereps)

    output:
    tuple val(readmode), path("all.dd.${readmode}.RDS"), emit: inferred

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def bandsize = params.platform == 'pacbio' ? ', BAND_SIZE=32' : ''
    def dadaOpt = params.dada_opts ? "${params.dada_opts}" : "NA"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    dadaOpt <- "${dadaOpt}"

    if (!is.na(dadaOpt)) {
      setDadaOpt(${dadaOpt})
      cat("dada Options:\\n",${dadaOpt},"\\n")
    }
    set.seed(100)

    cat("Processing all samples\\n")

    err <- readRDS("${err}")

    #Variable selection from CLI input flag --pool
    pool <- "${params.pool}"

    # 'pool' is a weird flag, either 'pseudo' (string), or T/F (bool)
    if(pool != "pseudo"){
      pool <- as.logical(pool)
    }
    
    # File parsing (these come from the process input channel)
    derep_files <- list.files('.', pattern=paste0("${readmode}",".derep.RDS"), full.names = TRUE)

    dereps <- lapply(derep_files, readRDS)

    # note this is a bit of a hack, but we want the file name 
    # included with the name of the derep object. This makes
    # sure these are in sync if needed later
    names(dereps) <- sapply(dereps, function(x) { x\$file })

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
