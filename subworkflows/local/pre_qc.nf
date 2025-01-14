include { FASTQC                 } from '../../modules/nf-core/fastqc/main'
include { PLOT_QUALITY_PROFILE   } from '../../modules/local/plotqualityprofile'
include { VSEARCH_EESTATS        } from '../../modules/local/vsearch_eestats'
include { VSEARCH_OVERLAP        } from '../../modules/local/vsearchoverlap'
include { MERGE_OVERLAP_CHECK    } from '../../modules/local/mergeoverlapcheck'
include { OVERLAP_HEATMAP        } from '../../modules/local/overlapheatmap'

workflow PRE_QC {

    take:
    ch_samplesheet
    skip_ee
    skip_merging
    skip_dadaQC

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    FASTQC (
        ch_samplesheet
    )
    
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    if (!skip_merging) {
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

    if (!skip_ee) {
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

