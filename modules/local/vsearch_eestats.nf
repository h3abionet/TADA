// TODO: at the moment I'm checking on the feasibility of using 
//       these for additional QC plots, esp on cum expected errors;
//       at the moment they are a bit of a data dump

process VSEARCH_EESTATS {
    tag "$meta.id"
    label 'process_low'

    container "quay.io/biocontainers/vsearch:2.27.0--h6a68c12_1"

    input:
    tuple val(meta), path(reads)

    output:
    path("${meta.id}.{R1,R2}.stats"), emit: stats
    path("${meta.id}.{R1,R2}.eestats"), emit: eestats
    path("${meta.id}.{R1,R2}.eestats2"), emit: eestats2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        vsearch --fastq_stats ${reads[0]} \
            --log ${meta.id}.R1.stats

        vsearch --fastq_eestats ${reads[0]} \
            --output ${meta.id}.R1.eestats

        vsearch --fastq_eestats2 ${reads[0]} \
            --ee_cutoffs "0.5,1.0,2.0,4.0,8.0" \
            --output ${meta.id}.R1.eestats2
        """
    } else {
        """
        vsearch --fastq_stats ${reads[0]} \
            --log ${meta.id}.R1.stats

        vsearch --fastq_stats ${reads[1]} \
            --log ${meta.id}.R2.stats            

        vsearch --fastq_eestats ${reads[0]} \
            --output ${meta.id}.R1.eestats

        vsearch --fastq_eestats ${reads[1]} \
            --output ${meta.id}.R2.eestats

        vsearch --fastq_eestats2 ${reads[0]}  \
            --ee_cutoffs "0.5,1.0,2.0,4.0,8.0" \
            --output ${meta.id}.R1.eestats2

        vsearch --fastq_eestats2 ${reads[1]} \
            --ee_cutoffs "0.5,1.0,2.0,4.0,8.0" \
            --output ${meta.id}.R2.eestats2
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.PDF
    """
}
