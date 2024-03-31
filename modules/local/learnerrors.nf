process LEARNERRORS {
    tag "$readmode"
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    tuple val(readmode), file(reads)

    output:
    // TODO nf-core: List additional required output channels/values here
    // path "versions.yml"           , emit: versions
    tuple val(readmode), file("errors.${readmode}.RDS"), emit: error_models
    // path("errors.R[12].RDS"), emit: errorModelsPerSample
    path("${readmode}.err.pdf"), emit: pdf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    def derepreads = 100000
    // TODO: not even sure the below will work w/ nf-core style parameters, 
    // may need to go to a string or args-type modifications
    def dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
    // TODO: will need to if-else this when implementing alternative error models
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    dadaOpt <- ${dadaOpt}

    if (!is.na(dadaOpt)) {
        setDadaOpt(dadaOpt)
        cat("dada Options:\n",dadaOpt,"\n")
    }

    # File parsing (these come from the process input channel)
    filts <- list.files('.', pattern=paste0("${readmode}",".filtered.fastq.gz"), full.names = TRUE)
    sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
    set.seed(100)

    # Learn forward error rates
    err <- learnErrors(filts, multithread=${task.cpus}, verbose=1)

    # This is a rough correction for NovaSeq binning issues
    # See https://github.com/h3abionet/TADA/issues/31, we'll likely
    # add alternatives here soon

    if (as.logical("${params.quality_binning}") == TRUE ) {
        # TODO: add alternatives for binning
        print("Running binning correction")
        errs <- t(apply(getErrors(err), 1, function(x) { x[x < x[40]] = x[40]; return(x)} ))
        err\$err_out <- errs
    }

    pdf(paste0("${readmode}",".err.pdf"))
    plotErrors(err, nominalQ=TRUE)
    dev.off()

    saveRDS(err, paste0("errors.","${readmode}",".RDS"))
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    # TODO: make a proper stub
    """
}
