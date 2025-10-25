process ILLUMINA_CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/cutadapt:5.0--py39hbcbf7aa_0'

    input:
    tuple val(meta), path(reads)
    val(for_primer)
    val(rev_primer)
    val(for_primer_rc)
    val(rev_primer_rc)

    output:
    tuple val(meta), path("${meta.id}.R[12].filtered.fastq.gz"), optional: true, emit: trimmed
    path("${meta.id}.cutadapt.out"), emit: trimmed_report // to merging data
    path("${meta.id}.cutadapt.json"), emit: cutadapt_json  // to MultiQC

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    trunc_for = params.trunc_for >=0 ? "-l ${params.trunc_for}" : ""
    trunc_rev = params.trunc_rev >=0 ? "-L ${params.trunc_rev}" : ""
    maxN = params.maxN >=0 ? "--max-n ${params.maxN}" : ""
    maxEE = "--max-ee ${[params.maxEE_for,params.maxEE_rev].max()}"
    min_len = params.min_read_len ? "-m ${params.min_read_len}" : "-m 50" 
    max_len = params.max_read_len != "Inf" ? "-M ${params.max_read_len}" : ""
    outr2 = meta.single_end ? '' : "-p ${meta.id}.R2.filtered.fastq.gz"
    p2 = meta.single_end ? '' : "-G ${rev_primer} -A ${for_primer_rc}"
    polyG = params.illumina_twocolor ? "--nextseq-trim=2" : ""

    """
    cutadapt \\
        --report=minimal \\
        --json=${meta.id}.cutadapt.json \\
        -g ${for_primer} -a ${rev_primer_rc} ${p2} \\
        --cores ${task.cpus} \\
        -n 2 ${maxEE} ${min_len} ${max_len} ${maxN} ${polyG} ${trunc_for} ${trunc_rev} \\
        -o ${meta.id}.R1.filtered.fastq.gz ${outr2} \\
        ${reads} > ${meta.id}.cutadapt.out

    # is the FASTQ file empty?
    if [ -n "\$(gunzip <${meta.id}.R1.filtered.fastq.gz | head -c 1 | tr '\\0\\n' __)" ]; then
        echo "Sequences present"
    else
        rm ${meta.id}.R[12].filtered.fastq.gz
    fi
    """

    // stub:
    // def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // """
    // """
}
