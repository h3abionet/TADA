process DADA2_LEARN_ERRORS {
    tag "${params.platform} ${readmode}"
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:docker-DADA-1.36"

    input:
    tuple val(readmode), path(reads)

    output:
    tuple val(readmode), file("errors.${readmode}.RDS"), emit: error_models
    path("${readmode}.err.pdf"), emit: pdf

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    def derepreads = 100000
    // TODO: 
    def dadaOpt = !params.dada_opts.isEmpty() ? "'${params.dada_opts.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
    if (params.platform == 'pacbio') {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        dadaOpt <- ${dadaOpt}

        if (!is.na(dadaOpt)) {
            setDadaOpt(dadaOpt)
            cat("dada Options:\n",dadaOpt,"\n")
        }
        
        # File parsing
        filts <- list.files('.', pattern="filtered.fastq.gz", full.names = TRUE)
        sample.namesF <- sapply(strsplit(basename(filts), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
        set.seed(100)

        dereps <- derepFastq(filts, n=${derepreads}, verbose=TRUE)

        # Learn forward error rates
        errs <- learnErrors(dereps, 
            nbases = 1e8,
            errorEstimationFunction=PacBioErrfun, 
            BAND_SIZE=32, 
            multithread=${task.cpus},
            verbose=TRUE)

        pdf(paste0("${readmode}",".err.pdf"))
        plotErrors(errs, nominalQ=TRUE)
        dev.off()

        saveRDS(errs, paste0("errors.","${readmode}",".RDS"))
        """
    } else if (params.platform == 'illumina') {
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
    }
}
