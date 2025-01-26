process VARIABLEFILTER {
    tag "$meta.id"
    label 'process_medium'
    container ""

    input:
    tuple val(meta), file(reads), file(trimming) from itsStep3.join(itsStep3Trimming)

    output:
    tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true, emit: filteredReadsR1
    tuple val(meta), file("${meta.id}.R2.filtered.fastq.gz") optional true, emit: filteredReadsR2
    tuple val(meta), file("${meta.id}.R[12].filtered.fastq.gz") optional true, emit: reads
    file "*.trimmed.txt", emit: read_tracking

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(ShortRead))
    suppressPackageStartupMessages(library(Biostrings))

    out <- filterAndTrim(fwd = paste0("${meta.id}",".R1.cutadapt.fastq.gz"),
                        filt = paste0("${meta.id}", ".R1.filtered.fastq.gz"),
                        rev = if("${reads[1]}" == "null") NULL else paste0("${meta.id}",".R2.cutadapt.fastq.gz"),
                        filt.rev = if("${reads[1]}" == "null") NULL else paste0("${meta.id}", ".R2.filtered.fastq.gz"),
                        maxEE = if("${reads[1]}" == "null") ${params.maxEEFor} else c(${params.maxEEFor}, ${params.maxEERev}), 
                        truncQ = ${params.truncQ},
                        rm.phix = as.logical(${params.rmPhiX}),
                        maxLen = ${params.max_read_len},
                        minLen = ${params.min_read_len},
                        compress = TRUE,
                        verbose = TRUE,
                        multithread = ${task.cpus})
    #Change input read counts to actual raw read counts
    colnames(out) <- c('cutadapt', 'filtered')
    write.csv(out3, paste0("${meta.id}", ".trimmed.txt"))
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${meta.id}.trimmed.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        variablefilter: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
