process VARIABLE_TRIM {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'

    input:
    tuple val(meta), path(reads)
    tuple val(for_primer), val(rev_primer)
    tuple val(for_primer_rc), val(rev_primer_rc)

    output:
    tuple val(meta), file("${meta.id}.R[12].cutadapt.fastq.gz") optional true, emit: trimmed_reads
    file("*.cutadapt.out") into cutadaptToMultiQC
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    outr2 = meta.single_end ? '' : "-p ${meta.id}.R2.cutadapt.fastq.gz"
    p2 = meta.single_end ? '' : "-G ${rev_primer} -A ${rev_primer_rc}"
    """
    cutadapt -g ${for_primer} -a ${for_primer_rc} ${p2} \\
        --cores ${task.cpus} \\
        --max-N ${params.maxN} \\
        -n 2 \\
        -o ${meta.id}.R1.cutadapt.fastq.gz ${outr2} \\
        ${reads} > ${meta.id}.cutadapt.out
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > "${meta.id}.R1.cutadapt.fastq.gz"
    touch ${meta.id}.cutadapt.out
    """
}
