process VSEARCH_OVERLAP {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/biocontainers/vsearch:2.27.0--h6a68c12_1"

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.id}.lengthstats.txt"), emit: merged_stats
    path("${meta.id}.log"),             emit: merged_log
    path "versions.yml",                emit: versions

    when:
    params.check_merging || task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    vsearch --fastq_mergepairs \\
        "${reads[0]}" \\
        --reverse "${reads[1]}" \\
        --fastqout "${meta.id}.merged.fastq" \\
        --threads ${task.cpus} \\
        --fastq_minovlen 5 \\
        --fastq_allowmergestagger \\
        --log "${meta.id}.log"

    awk 'NR%4 == 2 {lengths[length(\$0)]++} END {for (l in lengths) {print l, "\\t", lengths[l]}}' \\
        "${meta.id}.merged.fastq" \\
        > "${meta.id}.lengthstats.txt"

    rm "${meta.id}.merged.fastq"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearchoverlap: \$(vsearch --version |& sed -r '1!d ; s/vsearch\\s(\\S+),\\s+.*/\\1/')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch "${meta.id}.log"
    echo "" | gzip > ${meta.id}.merged.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearchoverlap: \$(vsearch --version |& sed -r '1!d ; s/vsearch\\s(\\S+),\\s+.*/\\1/')
    END_VERSIONS
    """
}
