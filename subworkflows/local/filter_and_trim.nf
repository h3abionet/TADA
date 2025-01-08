include { ILLUMINA_DADA2_FILTER_AND_TRIM   } from '../../modules/local/illumina_filterandtrim'
include { PACBIO_DADA2_FILTER_AND_TRIM     } from '../../modules/local/pacbio_filterandtrim'
include { PACBIO_CUTADAPT                  } from '../../modules/local/pacbio_cutadapt'
include { MERGE_TRIM_TABLES                } from '../../modules/local/mergetrimtables'

workflow FILTER_AND_TRIM {

    take:
    input //channel: [val(meta), path(reads)]

    main:
    // Three options for Illumina data:
    //       DADA2 trimming and filtering (PE and SE currently) - implemented
    //       cutadapt-based (primers + Ns) + vsearch (EE) - NYI
    //       Hybrid (variable length) - NYI
    // Two options for PacBio:
    //       DADA2 filtering (filter) - NYI
    //       cutadapt (trim) + vsearch - NYI

    ch_reports = Channel.empty()
    ch_trimmed = Channel.empty()
    ch_trimmed_R1 = Channel.empty()
    ch_trimmed_R2 = Channel.empty()

    if (params.platform == "pacbio") {

        // TODO: this could be modified/split into a `cutadapt`-only step; there
        // are additional filters for max EE and max N in cutadapt
        PACBIO_CUTADAPT(
            input
        )

        // TODO: should be summarized as well, go to MultiQC
        // ch_reports = PACBIO_CUTADAPT_FILTER_AND_TRIM.out.cutadapt_report.collect()
        
        // TODO: this could be modified/split into a `DADA2`-only step
        PACBIO_DADA2_FILTER_AND_TRIM(
            PACBIO_CUTADAPT.out.cutadapt_trimmed
        )
        ch_trimmed = PACBIO_DADA2_FILTER_AND_TRIM.out.trimmed
        ch_trimmed_R1 = PACBIO_DADA2_FILTER_AND_TRIM.out.trimmed
        ch_reports = PACBIO_DADA2_FILTER_AND_TRIM.out.trimmed_report.collect()
    } else {
        // this handles both paired and single-end data
        ILLUMINA_DADA2_FILTER_AND_TRIM(
            input
        )
        ch_trimmed = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed
        ch_reports = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_report.collect()
        ch_trimmed_R1 = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_R1
        ch_trimmed_R2 = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_R2
    }

    // TODO: add variable-length and PacBio
    MERGE_TRIM_TABLES(
        ch_reports
    )

    ch_trimmed_infer = ch_trimmed_R1
            .map { [ 'R1', it[1]] }
            .concat(ch_trimmed_R2.map {['R2', it[1]] } )
            .groupTuple(sort: true)
    // ch_trimmed_infer.dump(tag: "infer:", pretty: true)
    
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
    trimmed = ch_trimmed
    trimmed_report = MERGE_TRIM_TABLES.out.trimmed_report // channel: [ RDS ]
    trimmed_infer = ch_trimmed_infer
}
