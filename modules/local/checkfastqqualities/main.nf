process CHECK_FASTQ_QUALITIES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)

    output:
    // tuple val(meta), path("*.pdf"), emit: qc
    tuple val(meta), path("quality_bins.RDS"), emit: quality_bins_rds
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(ShortRead))
    suppressPackageStartupMessages(library(tidyverse))

    # Read just R1
    l <- FastqSampler(${reads[0]}, 
        n=1e6,
        readerBlockSize=1e4,

    fq <- yield(l)
    qual_matrix <- as(quality(fq), "matrix")
    qual_df <- as.data.frame(qual_matrix)
    qual_df$ReadID <- rownames(qual_df)
    qual_long <- qual_df %>%
      pivot_longer(
        cols = -ReadID,
        names_to = "BasePosition",
        values_to = "QualityScore"
      ) %>%
      mutate(BasePosition = as.integer(gsub("V", "", BasePosition)))

    # this assumes there are 10 or fewer bins
    stopifnot(length(binned_quals) <= 10)

    binned_quals <- factor(qual_long$QualityScore) %>% 
      levels() %>% 
      as.integer()

    saveRDS(binned_quals, "quality_bins.RDS")
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkfastqqualities: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
