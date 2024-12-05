include { ILLUMINA_FILTER_AND_TRIM   } from '../../modules/local/filterandtrim'
include { PACBIO_FILTER_AND_TRIM     } from '../../modules/local/filterandtrim'
include { MERGE_TRIM_TABLES          } from '../../modules/local/mergetrimtables'

workflow FILTER_AND_TRIM {

    take:
    input //channel: [val(meta), path(reads)]

    main:
    // Three options for Illumina data:
    //       DADA2 trimming and filtering (PE and SE currently) - implemented
    //       cutadapt-based (primers + Ns) + vsearch (EE) - NYI
    //       Hybrid (variable length) - NYI
    // Two options for PacBio:
    //       cutadapt (trim) + DADA2 filtering (filter) - NYI
    //       cutadapt (trim) + vsearch - NYI

    // DADA2

    // if (platform == 'pacbio' || platform == 'pacbio-kinnex') {
    // process PacBioTrim {
    //     tag { "PacBioTrim_${meta.id}" }
    //     publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

    //     input:
    //     // TODO: Note the channel name here should probably be changed
    //     tuple val(meta), file(reads) from dada2ReadPairs

    //     output:
    //     // tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1
    //     tuple val(meta), file("${meta.id}.noprimer.fastq.gz") optional true into filteredReadsToDADA2
    //     file("*.cutadapt.out") into cutadaptToMultiQC
    //     file("${meta.id}.untrimmed.fastq.gz")

    //     when:
    //     !(params.precheck)

    //     script:
    //     strictness = params.pacbio_strict_match ? '-g' : '-a'
    //     """
    //     # Logic: we should trim out the HiFi reads and require *both* primers be present (-g).
    //     # This should also reorient the sequence to match the primers (--rc).
    //     # Keep anything longer than 50bp, and allow users to filter their data by length later
    //     revprimer_rc=\$( echo -n ${params.revprimer} | tr "[ATGCUNYRSWKMBDHV]" "[TACGANRYSWMKVHDB]" | rev )

    //     cutadapt --rc \\
    //         ${strictness} "${params.fwdprimer}...\${revprimer_rc}" \\
    //         -m 50 \\
    //         -j ${task.cpus} \\
    //         --untrimmed-output "${meta.id}.untrimmed.fastq.gz" \\
    //         -o "${meta.id}.noprimer.fastq.gz" \\
    //         ${reads} > "${meta.id}.noprimer.cutadapt.out"
    //     """
    // }

    // process PacBioFilter {
    //     tag { "PacBioFilter_${meta.id}" }
    //     publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

    //     input:
    //     // TODO: Note the channel name here should probably be changed
    //     tuple val(meta), file(reads) from filteredReadsToDADA2

    //     output:
    //     // tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1
    //     tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1,readsToFastQC,readsToPerSample
    //     file "*.trimmed.txt" into trimTracking

    //     when:
    //     !(params.precheck)

    //     script:
    //     template "PacBioFilterAndTrim.R"
        
    // }

    // filteredReadsR2 = Channel.empty()

    // } else if (platform == 'illumina' && ( params.amplicon == 'ITS' || params.amplicon == 'variable'))  {

        // // this path is only needed when using variable length sequences
        // process ITSFilterAndTrimStep1 {
        //     tag { "ITS_Step1_${meta.id}" }

        //     input:
        //     tuple val(meta), file(reads) from dada2ReadPairs

        //     output:
        //     tuple val(meta), file("${meta.id}.R[12].noN.fastq.gz") optional true into itsStep2
        //     tuple val(meta), file("${meta.id}.out.RDS") into itsStep3Trimming  // needed for join() later
        //     file('forward_rc') into forwardP
        //     // TODO make this optional if data are SE
        //     file('reverse_rc') into reverseP

        //     when:
        //     !(params.precheck)

        //     script:
        //     template "ITSFilterAndTrimStep1.R"
        // }
        
        // process ITSFilterAndTrimStep2 {
        //     tag { "ITS_Step2_${meta.id}" }
        //     publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

        //     input:
        //     tuple(meta), file(reads) from itsStep2
        //     file(forP) from forwardP
        //     file(revP) from reverseP
            
        //     output:
        //     tuple val(meta), file("${meta.id}.R[12].cutadapt.fastq.gz") optional true into itsStep3
        //     file("*.cutadapt.out") into cutadaptToMultiQC

        //     when:
        //     !(params.precheck)

        //     script:
        //     outr2 = meta.single_end ? '' : "-p ${meta.id}.R2.cutadapt.fastq.gz"
        //     p2 = meta.single_end ? '' : "-G ${params.revprimer} -A \$REV_PRIMER"
        //     """
        //     FWD_PRIMER=\$(<forward_rc)
        //     REV_PRIMER=\$(<reverse_rc)
            
        //     cutadapt -g ${params.fwdprimer} -a \$FWD_PRIMER ${p2} \\
        //         --cores ${task.cpus} \\
        //         -n 2 \\
        //         -o ${meta.id}.R1.cutadapt.fastq.gz ${outr2} \\
        //         ${reads} > ${meta.id}.cutadapt.out
        //     """
        // }

        // process ITSFilterAndTrimStep3 {
        //     tag { "ITS_Step3_${meta.id}" }
        //     publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

        //     input:
        //     tuple val(meta), file(reads), file(trimming) from itsStep3.join(itsStep3Trimming)

        //     output:
        //     tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1
        //     tuple val(meta), file("${meta.id}.R2.filtered.fastq.gz") optional true into filteredReadsR2
        //     tuple val(meta), file("${meta.id}.R[12].filtered.fastq.gz") optional true into readsToFastQC,readsToPerSample
        //     file "*.trimmed.txt" into trimTracking

        //     when:
        //     !(params.precheck)

        //     script:
        //     template "ITSFilterAndTrimStep3.R"
        // }
    // else if (platform == 'illumina' && params.amplicon == '16S'){
    ILLUMINA_FILTER_AND_TRIM(
        input
    )

    // } else {
    //     // We need to shut this down!
    //     Channel.empty().into {cutadaptToMultiQC;filteredReads;filteredReadsforQC}
    // }

    ch_reports = ILLUMINA_FILTER_AND_TRIM.out.trimmed_report.collect()

    // TODO: add variable-length and PacBio
    MERGE_TRIM_TABLES(
        ch_reports
    )

    // Channel setup

    // We need to group data depending on which downstream steps are needed.  There
    // are two combinations possible

    // 1. The immediate downstream QC steps can use the meta info and the read pairs.
    //    Instead of doing handstands reusing the two channels above, we emit channels 
    //    with the reads paired if needed.

    // 2. LearnErrors and the pooled denoising branch requires all R1 and all R2, but 
    //    the two groups can be processed in parallel.  So we set up the channels with 
    //    this in mind. No sample ID info is really needed.
    // ch_trimmed_infer = FILTERANDTRIM.out.trimmed_R1
    //         .map { [ 'R1', it[1]] }
    //         .concat(FILTERANDTRIM.out.trimmed_R2.map {['R2', it[1]] } )
    //         .groupTuple(sort: true)
    emit:
    trimmed = ILLUMINA_FILTER_AND_TRIM.out.trimmed
    trimmed_report = MERGE_TRIM_TABLES.out.trimmed_report // channel: [ RDS ]
    trimmed_infer = ILLUMINA_FILTER_AND_TRIM.out.trimmed_R1
            .map { [ 'R1', it[1]] }
            .concat(ILLUMINA_FILTER_AND_TRIM.out.trimmed_R2.map {['R2', it[1]] } )
            .groupTuple(sort: true)
    // versions = ch_versions                     // channel: [ versions.yml ]
}

