process QIIME2_TREE {

    container "quay.io/qiime2/core:2021.4"

    input:
    path(unrooted)
    path(rooted)

    output:
    path("unrooted_tree.qza")
    path("rooted_tree.qza")

    when:
    task.ext.when == null || task.ext.when 

    script:
    def args = task.ext.args ?: ''
    """
    qiime tools import \
        --input-path ${unrooted} \
        --output-path unrooted_tree.qza \
        --type 'Phylogeny[Unrooted]'

    qiime tools import \
        --input-path ${rooted} \
        --output-path rooted_tree.qza \
        --type 'Phylogeny[Rooted]'
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
