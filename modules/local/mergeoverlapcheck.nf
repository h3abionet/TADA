process MERGE_OVERLAP_CHECK {
    tag 'Overlap Stats'
    label 'process_single'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(logs)

    output:
    path("all_mergedstats.tsv"), emit: overlap_check_stats
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    MergeCheck_Table.R
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}