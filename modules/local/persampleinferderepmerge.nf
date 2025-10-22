process PER_SAMPLE_INFER {
    tag "$meta.id"
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:docker-DADA-1.36"

    input:
    tuple val(meta), path(reads)
    path(errs)
    // optional inputs
    path(fp, stageAs: "priors_R1")
    path(rp, stageAs: "priors_R2")

    output:
    path("${meta.id}.{R1,merged}.RDS"), emit: combinedReads
    tuple val(meta), path("${meta.id}.dd.R{1,2}.RDS"), emit: dds
    val(readmode), emit: readmode

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bandsize = params.platform == 'pacbio' ? ', BAND_SIZE=32' : ''
    def dadaOpt = !params.dada_opts.isEmpty() ? "'${params.dada_opts.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
    readmode = errs.size() == 2 ? 'merged' : 'R1'
    run_fpriors = params.for_priors ? "TRUE" : "FALSE"
    run_rpriors = params.rev_priors ? "TRUE" : "FALSE"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(Biostrings))
    suppressPackageStartupMessages(library(tidyverse))

    set.seed(100)

    getPriors <- function(x) {
      priors <- readDNAStringSet(x) |> as.vector() |> unname()
      return(priors)
    }

    dadaOpt <- ${dadaOpt}

    if (!is.na(dadaOpt)) {
      setDadaOpt(dadaOpt)
      cat("dada Options:\\n",dadaOpt,"\\n")
    }

    cat("Processing:", "${meta.id}", "\\n")

    errF <- readRDS("errors.R1.RDS")
    derepF <- derepFastq("${reads[0]}", n=100000)

    # TODO: there is probably a better way of doing this 
    # when using optparse
    paramsF <- list(
        derep=derepF, 
        err=errF,
        multithread=${task.cpus},
        pool=FALSE ${bandsize}
    )

    if (as.logical("${run_fpriors}")) {
      paramsF\$priors <- getPriors("${fp}")
    }

    ddF <- do.call(dada, paramsF)
    saveRDS(ddF, "${meta.id}.dd.R1.RDS")

    if (file.exists("errors.R2.RDS")) {
        errR <- readRDS("errors.R2.RDS")
        derepR <- derepFastq("${reads[1]}", n=100000)
        paramsR <- list(
            derep=derepR, 
            err=errR, 
            multithread=${task.cpus}, 
            pool=FALSE ${bandsize}
        )

        if (as.logical("${run_rpriors}")) {
            paramsR\$priors <- getPriors("${rp}")
        }

        message("DADA2 params, R2:", paramsR, "\\n")
        ddR <- do.call(dada, paramsR)
        saveRDS(ddR, "${meta.id}.dd.R2.RDS")

        merger <- mergePairs(ddF, derepF, ddR, derepR,
            returnRejects = TRUE,
            minOverlap = ${params.min_overlap},
            maxMismatch = ${params.max_mismatch},
            trimOverhang = as.logical("${params.overhang_trim}"),
            justConcatenate=as.logical("${params.just_concatenate}")
        )

        saveRDS(merger, paste("${meta.id}.merged.RDS", sep="."))
    } else {
        # yes this is a little silly (it's the same as the dd.R1.RDS above).
        # But it does make the logical flow through this channel easier
        saveRDS(ddF, paste("${meta.id}.R1.RDS", sep="."))
    }
    """

    // stub:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // """
    // """
}
