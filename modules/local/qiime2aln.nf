process QIIME2_ALIGNMENT {

    container "quay.io/qiime2/core:2021.4"

    input:
    path(aln)

    output:
    path("aligned_asvs.qza")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    qiime tools import \
        --input-path ${aln} \
        --output-path aligned_asvs.qza \
        --type 'FeatureData[AlignedSequence]'
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
