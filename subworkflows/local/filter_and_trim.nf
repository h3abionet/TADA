include { ILLUMINA_FILTER_AND_TRIM   } from '../../modules/local/filterandtrim'
// include { PACBIO_FILTER_AND_TRIM     } from '../../modules/local/filterandtrim'
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

    ILLUMINA_FILTER_AND_TRIM(
        input
    )

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
}

