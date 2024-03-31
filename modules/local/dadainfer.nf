process DADAINFER {
    tag "$readmode"    
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    tuple val(readmode), file(err), file(reads)

    output:
    // TODO nf-core: List additional required output channels/values here
    // path "versions.yml"           , emit: versions
    path("all.dd.${readmode}.RDS"), emit: inferred

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    dadaOpt <- ${dadaOpt}

    if (!is.na(dadaOpt)) {
      setDadaOpt(dadaOpt)
      cat("dada Options:\\n",dadaOpt,"\\n")
    }

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
    dds <- dada(dereps, err=err, multithread=${task.cpus}, pool=pool)

    saveRDS(dds, "all.dd.${readmode}.RDS")
    """

    stub:
    def args = task.ext.args ?: ''
    
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    # add some real stuff here
    touch all.dd.${readmode}.RDS
    """
}
