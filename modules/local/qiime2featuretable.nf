process QIIME2_FEATURETABLE {

    container "quay.io/qiime2/core:2021.4"
    
    input:
    path(seqtab)

    output:
    path("seqtab_final.qza")

    when:
    task.ext.when == null || task.ext.when 

    script:
    def args = task.ext.args ?: ''
    """
    biom convert -i ${seqtab} \
        -o seqtab-biom-table.biom \
        --table-type="OTU table" \
        --to-hdf5

    qiime tools import \
        --input-path seqtab-biom-table.biom \
        --input-format BIOMV210Format \
        --output-path seqtab_final.qza \
        --type 'FeatureTable[Frequency]'
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    """
}
