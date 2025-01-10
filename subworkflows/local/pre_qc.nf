// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { PLOT_QUALITY_PROFILE   } from '../../modules/local/plotqualityprofile'
include { VSEARCH_EESTATS        } from '../../modules/local/vsearch_eestats'
include { VSEARCH_OVERLAP        } from '../../modules/local/vsearchoverlap'
include { MERGE_OVERLAP_CHECK    } from '../../modules/local/mergeoverlapcheck'
include { OVERLAP_HEATMAP        } from '../../modules/local/overlapheatmap'

workflow PRE_QC {

    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()

    VSEARCH_OVERLAP(
        ch_samplesheet
    )

    ch_versions = ch_versions.mix(VSEARCH_OVERLAP.out.versions.first())

    VSEARCH_EESTATS (
        ch_samplesheet
    )

    ch_versions = ch_versions.mix(VSEARCH_EESTATS.out.versions.first())

    MERGE_OVERLAP_CHECK(
        VSEARCH_OVERLAP.out.merged_log.collect()
    )

    OVERLAP_HEATMAP(
        VSEARCH_OVERLAP.out.merged_stats.collect()
    )

    PLOT_QUALITY_PROFILE (
        ch_samplesheet
    )

    ch_versions = ch_versions.mix(PLOT_QUALITY_PROFILE.out.versions.first())
    
    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}
