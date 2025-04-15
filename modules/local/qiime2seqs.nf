process QIIME2_SEQUENCE {

    container "quay.io/qiime2/core:2021.4"

    input:
    path(seqs)

    output:
    path("asv_sequences.qza")

    when:
    task.ext.when == null || task.ext.when 

    script:
    def args = task.ext.args ?: ''
    """
    qiime tools import \
        --input-path ${seqs} \
        --output-path asv_sequences.qza \
        --type 'FeatureData[Sequence]'
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
