process OVERLAP_HEATMAP {
    tag "Overlap Heatmap"
    label 'process_single'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(merged_tables)

    output:
    path("MergedCheck_heatmap.pdf"), emit: overlap_check_pdf
    path("MergedCheck_heatmap.RDS"), emit: overlap_check_rds
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    MergeCheck_Plot.R
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
