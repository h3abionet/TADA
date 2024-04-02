// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { ILLUMINA_FILTER_AND_TRIM   } from '../../modules/local/filterandtrim'
include { MERGE_TRIM_TABLES          } from '../../modules/local/mergetrimtables'

workflow FILTER_AND_TRIM {

    take:
    input //channel: [val(meta), path(reads)

    main:
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
    // versions = ch_versions                     // channel: [ versions.yml ]
}

