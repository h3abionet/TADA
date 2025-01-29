/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/nf-core/multiqc/main'

include { PRE_QC                 } from '../subworkflows/local/pre_qc'
include { FILTER_AND_TRIM        } from '../subworkflows/local/filter_and_trim'
include { DADA2_DENOISE          } from '../subworkflows/local/dada2_denoise'
include { TAXONOMY               } from '../subworkflows/local/taxonomy'
include { PHYLOGENY              } from '../subworkflows/local/phylogeny'
include { QUALITY_CONTROL        } from '../subworkflows/local/qualitycontrol'
include { GENERATE_OUTPUT        } from '../subworkflows/local/generate_output'

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_tada_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TADA {

    take:
    ch_samplesheet // channel: [val(meta), path(reads)]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_readtracking = Channel.empty()

    // Make sure this is used throughout the workflow instead of params.platform
    def platform = params.platform.toLowerCase()

    if (!(["illumina","pacbio"].contains(platform))) {
        exit 1, "Only supported platforms (--platform argument) are currently 'pacbio' or 'illumina'"
    }

    // TODO: implement seqtable input?
    // // ${deity} there has to be a better way to check this!
    // if ( (params.seqTables && params.reads) || 
    //      (params.input  && params.reads) || 
    //      (params.seqTables && params.input)) {
    //     exit 1, "Only one of --reads, --input, or --seqTables is allowed!"
    // }

    // TODO: rethink this step and what the expected inputs mean
    // if ( params.trimming && (!(params.fwd_adaptor) || !(params.rev_adapter))) ) {
    //     log.info "Both --fwdprimer and --revprimer should be set unless skipping all trimming steps.\n" 
    //     log.info "[These options will become requirements in future releases]"
    // } else { 
    //     // this is a check if the primers are supplied but trimming is also set, which is almost certainly a mistake
    //     if ( params.trim_for != 0 || params.trim_rev != 0 ) {
    //         log.info "trim_for and/or trim_rev are set along with one or more primers (fwd_adapter/rev_adapter).\n" 
    //         log.info "This will trim additional bases *in addition to* the supplied primers!"
    //     }
    // }

    // if (params.aligner == 'infernal' && params.infernalCM == false){
    //     exit 1, "Must set covariance model using --infernalCM when using Infernal"
    // }

    if (!(['simple','md5'].contains(params.id_type))) {
        exit 1, "--id_type can currently only be set to 'simple' or 'md5', got ${params.id_type}"
    }

    PRE_QC(
        ch_samplesheet,
        params.skip_FASTQC,
        params.skip_dadaQC,
        params.skip_merging_check,
        params.skip_ee_check
    )

    ch_multiqc_files = ch_multiqc_files.mix(PRE_QC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(PRE_QC.out.versions)
    
    // ch_multiqc_files = ch_multiqc_files.mix(PLOTQUALITYPROFILE.out.zip.collect{it[1]})

    FILTER_AND_TRIM (
        ch_samplesheet,
        params.skip_trimming
    )
    ch_multiqc_files = ch_multiqc_files.mix(FILTER_AND_TRIM.out.ch_multiqc_files)
    ch_readtracking = ch_readtracking.mix(FILTER_AND_TRIM.out.trimmed_report)

    // TODO: Input for these should be the trimmed reads from above, but
    // they may need to be mixed in different ways depending on the 
    // denoising workflow used
    
    // TODO: harmonize output when possible (FASTA + TSV)
    DADA2_DENOISE(
        FILTER_AND_TRIM.out.trimmed_infer
    )

    // TODO: split out chimera removal into a separate step
    // CHIMERA_REMOVAL()

    // TODO: mix these in a specific order? This would help when merging
    //       multiple tables from different sources
    ch_readtracking = ch_readtracking.mix(DADA2_DENOISE.out.inferred.collect())
    ch_readtracking = ch_readtracking.mix(DADA2_DENOISE.out.seqtable_renamed)
    ch_readtracking = ch_readtracking.mix(DADA2_DENOISE.out.merged_seqs)

    // Subworkflows-Taxonomic assignment (optional)
    ch_taxtab = Channel.empty()
    ch_boots =  Channel.empty()

    if (params.reference) {
        ref_file = file(params.reference, checkIfExists: true)
        species_file = params.species ? file(params.species, checkIfExists: true) : file("${projectDir}/assets/dummy_file")

        // TODO: add alternative callers:
        //      DECIPHER, USEARCH/VSEARCH,q2, BLAST
        // TODO: readmap -> FASTA?
        TAXONOMY(
            DADA2_DENOISE.out.readmap,
            ref_file,
            species_file
        )
        ch_taxtab = TAXONOMY.out.ch_taxtab
        ch_metrics = TAXONOMY.out.ch_metrics
    }
    
    PHYLOGENY(DADA2_DENOISE.out.nonchimeric_asvs)

    // Post-QC
    QUALITY_CONTROL(
        ch_readtracking,
        DADA2_DENOISE.out.merged_seqs,
        DADA2_DENOISE.out.filtered_seqtable
    )
    // READ_TRACKING(
    //     ch_readtracking.collect()
    // )

    // PLOT_MERGED_HEATMAP(
    //     DADA2_DENOISE.out.merged_seqs
    // )

    // PLOT_ASV_DIST(
    //     DADA2_DENOISE.out.filtered_seqtable
    // )

    GENERATE_OUTPUT(
        DADA2_DENOISE.out.seqtable_renamed,
        DADA2_DENOISE.out.seqtab2qiime,
        TAXONOMY.out.ch_taxtab_rds,
        TAXONOMY.out.ch_taxtab,
        DADA2_DENOISE.out.nonchimeric_asvs,
        PHYLOGENY.out.ch_alignment,
        PHYLOGENY.out.ch_unrooted_tree,
        PHYLOGENY.out.ch_rooted_tree
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
