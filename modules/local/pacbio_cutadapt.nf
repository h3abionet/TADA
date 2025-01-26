// TODO: this is currently a local module; we should try to set this up
// to use the standard nf-core module
process PACBIO_CUTADAPT {
    tag "${meta.id}"

    container 'quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'

    input:
    // TODO: Note the channel name here should probably be changed
    tuple val(meta), path(reads)

    output:
    tuple val(meta), file("${meta.id}.noprimer.fastq.gz"), optional: true, emit: cutadapt_trimmed
    // file("${meta.id}.cutadapt.out"), emit: cutadapt_report
    // file("${meta.id}.untrimmed.fastq.gz"), emit: cutadapt_untrimmed

    // when:
    // !(params.precheck)

    script:
    strictness = params.pacbio_strict_match ? '-g' : '-a'
    """
    # Logic: we should trim out the HiFi reads and require *both* primers be present (-g).
    # This should also reorient the sequence to match the primers (--rc).
    # Keep anything longer than 50bp, and allow users to filter their data by length later
    revprimer_rc=\$( echo -n ${params.rev_primer} | tr "[ATGCUNYRSWKMBDHV]" "[TACGANRYSWMKVHDB]" | rev )

    cutadapt --rc \\
        ${strictness} "${params.for_primer}...\${revprimer_rc}" \\
        -m 50 \\
        -j ${task.cpus} \\
        --untrimmed-output "${meta.id}.untrimmed.fastq.gz" \\
        -o "${meta.id}.noprimer.fastq.gz" \\
        ${reads} > "${meta.id}.cutadapt.out"
    """
}
