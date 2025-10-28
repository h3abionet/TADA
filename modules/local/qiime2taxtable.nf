process QIIME2_TAXTABLE {

    container "quay.io/qiime2/amplicon:2025.7"
    
    input:
    path(taxtab)

    output:
    path("taxtab.qza"), emit: taxtab_qza
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when 

    script:
    def args = task.ext.args ?: ''
    """
    tail -n +2 ${taxtab} | \
        perl -ne 's/\\"//g; @foo = split; print "\$foo[0]\\t".join(";", @foo[1..\$#foo])."\\n"' \
        > headerless.txt
    
    qiime tools import \
        --input-path headerless.txt \
        --input-format HeaderlessTSVTaxonomyFormat \
        --output-path taxtab.qza \
        --type 'FeatureData[Taxonomy]'

    qiime metadata tabulate \
      --m-input-file taxtab.qza \
      --o-visualization taxtab.qzv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch taxtab.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS    
    """
}
