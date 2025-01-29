include { FASTQC                 } from '../../modules/nf-core/fastqc/main'
include { PLOT_QUALITY_PROFILE   } from '../../modules/local/plotqualityprofile'
include { VSEARCH_EESTATS        } from '../../modules/local/vsearch_eestats'
include { VSEARCH_OVERLAP        } from '../../modules/local/vsearch_overlap'
include { MERGE_OVERLAP_CHECK    } from '../../modules/local/mergeoverlapcheck'
include { OVERLAP_HEATMAP        } from '../../modules/local/overlapheatmap'

workflow PRE_QC {

    take:
    ch_samplesheet
    skip_FASTQC
    skip_dadaQC
    skip_merging_check
    skip_ee_check

    main:
    ch_versions = Channel.empty()

    // TODO: this may need to be reimplemented if we add
    //       other outputs to MultiQC (e.g. custom ones)
    // ch_multiqc_files = Channel.empty()

    if (!skip_FASTQC) {
        FASTQC (
            ch_samplesheet
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }
    
    // TODO: this might be a step we run optionally after trimming (maybe a subworkflow?).
    //       This also doesn't need all data, maybe a random sampling would be best
    if (!skip_merging_check) {
        VSEARCH_OVERLAP(
            ch_samplesheet
        )

        ch_versions = ch_versions.mix(VSEARCH_OVERLAP.out.versions.first())

        MERGE_OVERLAP_CHECK(
            VSEARCH_OVERLAP.out.merged_log.collect()
        )

        OVERLAP_HEATMAP(
            VSEARCH_OVERLAP.out.merged_stats.collect()
        )        
    }

    if (!skip_ee_check) {
        VSEARCH_EESTATS (
            ch_samplesheet
        )

        ch_versions = ch_versions.mix(VSEARCH_EESTATS.out.versions.first())
    }
    
    if (!skip_dadaQC) {
        PLOT_QUALITY_PROFILE (
            ch_samplesheet
        )

        ch_versions = ch_versions.mix(PLOT_QUALITY_PROFILE.out.versions.first())
    }
    
    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
    zip = FASTQC.out.zip
}

