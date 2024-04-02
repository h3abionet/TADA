process MERGE_TRIM_TABLES {
    label 'process_low'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(trimData)

    output:
    path("all.trimmed.csv"), emit: trimmed_report

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
    """
    touch "all.trimmed.csv"
    """
}
