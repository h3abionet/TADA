// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { MMSEQS_CONVERTALIS } from '../../../modules/local/mmseqs/convertalis'
include { MMSEQS_CONVERTMSA  } from '../../../modules/local/mmseqs/convertmsa'
include { MMSEQS_MSA2PROFILE } from '../../../modules/local/mmseqs/msa2profile'
include { MMSEQS_SEARCH      } from '../../../modules/nf-core/mmseqs/search'
include { MMSEQS_CREATEDB    as MMSEQS_CREATEDB_PROFILE } from '../../../modules/nf-core/mmseqs/createdb'
include { MMSEQS_CREATEDB    as MMSEQS_CREATEDB_QUERY   } from '../../../modules/nf-core/mmseqs/createdb'
// include { MMSEQS_FILTER_DATA }

workflow MMSEQS_PROFILE_FILTER {

    take:
    // TODO nf-core: edit input (take) channels
    ch_seqtab  // channel: [ val(meta), [ bam ] ]
    ch_asvs    // channel: [ val(meta), [ bam ] ]
    ch_profile // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    MMSEQS_CONVERTMSA(
        ch_profile
    )

    MMSEQS_CREATEDB_PROFILE(
        MMSEQS_CONVERTMSA.out
    )

    MMSEQS_CREATEDB_QUERY(
        ch_asvs
    )

    MMSEQS_SEARCH(
        ch_query = MMSEQS_CREATEDB_QUERY.out,
        ch_db    = MMSEQS_CREATEDB_PROFILE.out,
    )

    MMSEQS_CONVERTALIS(
        MMSEQS_SEARCH.out
    )


    // TODO nf-core: substitute modules here for the modules of your subworkflow

    

    emit:
    // TODO nf-core: edit emitted channels
    // bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    // bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    // csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

