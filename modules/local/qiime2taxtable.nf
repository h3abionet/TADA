process QIIME2_TAXTABLE {

    container "quay.io/qiime2/core:2021.4"
    
    input:
    path(taxtab)

    output:
    path("taxtab_final.qza")

    when:
    task.ext.when == null || task.ext.when 

    script:
    def args = task.ext.args ?: ''
    """
    tail -n +2 ${taxtab} > headerless.txt
    qiime tools import \
        --input-path headerless.txt \
        --input-format HeaderlessTSVTaxonomyFormat \
        --output-path taxtab_final.qza \
        --type 'FeatureData[Taxonomy]'
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
