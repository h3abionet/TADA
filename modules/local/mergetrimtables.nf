process MERGETRIMTABLES {
    label 'process_low'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(trimData)

    output:
    path("all.trimmed.csv"), emit: trimmed_report
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    // tuple val(meta), path("*.bam"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    """
    #!/usr/bin/env Rscript

    # Main purpose of this script is to merge all trimming data into one table

    trimmedFiles <- list.files(path = '.', pattern = '*.trimmed.txt')
    sample.names <- sub('.trimmed.txt', '', trimmedFiles)
    trimmed <- do.call("rbind", lapply(trimmedFiles, function (x) as.data.frame(read.csv(x))))
    colnames(trimmed)[1] <- "Sequence"
    trimmed\$SampleID <- sample.names
    write.csv(trimmed, "all.trimmed.csv", row.names = FALSE)
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch "all.trimmed.csv"
    """
}
