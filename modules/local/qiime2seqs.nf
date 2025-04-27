process QIIME2_SEQUENCE {

    container "quay.io/qiime2/amplicon:2024.10"

    input:
    path(seqs)

    output:
    path("asv_sequences.qza"), emit: asvs_qza
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when 

    script:
    def args = task.ext.args ?: ''
    """
    qiime tools import \
        --input-path ${seqs} \
        --output-path asv_sequences.qza \
        --type 'FeatureData[Sequence]'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch asv_sequences.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS    
    """
}
