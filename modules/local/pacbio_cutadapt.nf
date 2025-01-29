process PACBIO_CUTADAPT {
    tag "${meta.id}"

    container 'quay.io/biocontainers/cutadapt:5.0--py39hbcbf7aa_0'

    input:
    tuple val(meta), path(reads)
    val(for_primer)
    val(rev_primer_rc)

    output:
    tuple val(meta), file("${meta.id}.filtered.fastq.gz"), optional: true, emit: cutadapt_trimmed
    file("${meta.id}.cutadapt.out"), emit: trimmed_report // to merging data
    file("${meta.id}.untrimmed.fastq.gz"), emit: cutadapt_untrimmed
    file("${meta.id}.cutadapt.json"), emit: cutadapt_json  // to MultiQC

    when:
    task.ext.when == null || task.ext.when

    script:
    strictness = params.pacbio_strict_match ? '-g' : '-a'
    maxN = params.maxN >=0 ? "--max-n ${params.maxN} " : ""
    maxEE = [params.max_ee_for,params.max_ee_rev].max() == 0 ? "--max-ee ${[params.max_ee_for,params.max_ee_rev].max()}" : ""
    min_len = params.min_read_len ? "-m ${params.min_read_len}" : "-m 50" 
    max_len = params.max_read_len != "Inf" ? "-M ${params.max_read_len}" : ""
    """
    cutadapt --rc \\
        --report=minimal \\
        ${strictness} "${for_primer}...${rev_primer_rc}" \\
        -j ${task.cpus} ${min_len} ${max_len} ${maxEE} ${max_N} \\
        --untrimmed-output "${meta.id}.untrimmed.fastq.gz" \\
        --json=${meta.id}.cutadapt.json \\
        -o "${meta.id}.filtered.fastq.gz" \\
        ${reads} > "${meta.id}.cutadapt.out"

    # is the FASTQ file empty?
    if [ -n "\$(gunzip <${meta.id}.filtered.fastq.gz | head -c 1 | tr '\\0\\n' __)" ]; then
        echo "Sequences present"
    else
        rm ${meta.id}.filtered.fastq.gz
    fi
    """
}
