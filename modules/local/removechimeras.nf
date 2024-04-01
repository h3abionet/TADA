process REMOVECHIMERAS {
    label 'process_medium'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(st)

    output:
    path("seqtab.nonchim.RDS"), emit: nonchim_seqtable

    when:
    task.ext.when == null || task.ext.when

    script:
    chimOpts = params.removeBimeraDenovo_options ? ", ${params.removeBimeraDenovo_options}" : ""
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    st.all <- readRDS("${st}")

    # Remove chimeras
    seqtab <- removeBimeraDenovo(
        st.all, 
        method="consensus", 
        multithread=${task.cpus}, 
        verbose=TRUE ${chimOpts} 
        )

    saveRDS(seqtab, "seqtab.nonchim.RDS")
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    """

    """
}
