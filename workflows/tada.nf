/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PLOT_QUALITY_PROFILE   } from '../modules/local/plotqualityprofile'
include { VSEARCH_EESTATS        } from '../modules/local/vsearch_eestats'

include { FILTER_AND_TRIM        } from '../subworkflows/local/filter_and_trim'
include { DADA2_DENOISE          } from '../subworkflows/local/dada2_denoise'

// TODO: may want to move into a subworkflow since we will likely implement a 
//       few additional methods (q2-feature-classifier, IDTAXA, etc)
include { ASSIGN_TAXA_SPECIES    } from '../modules/local/assigntaxaspecies'

// TODO: Move into phylogenetic subworkflow
include { DECIPHER               } from '../modules/local/decipher'
include { PHANGORN               } from '../modules/local/phangorn'
include { FASTTREE               } from '../modules/local/fasttree'
include { ROOT_TREE              } from '../modules/local/roottree'

// TODO: Move into subworkflow(s)
include { READ_TRACKING          } from '../modules/local/readtracking'
include { PLOT_MERGED_HEATMAP    } from '../modules/local/plotmerged'
include { PLOT_ASV_DIST          } from '../modules/local/plotasvlen'

// TODO: Move into subworkflow(s)
include { BIOM                   } from '../modules/local/biom'
include { SEQTABLE2TEXT          } from '../modules/local/seqtable2txt'
include { TAXTABLE2TEXT          } from '../modules/local/taxtable2txt'
include { QIIME2_FEATURETABLE    } from '../modules/local/qiime2featuretable'
include { QIIME2_TAXTABLE        } from '../modules/local/qiime2taxtable'
include { QIIME2_SEQUENCE        } from '../modules/local/qiime2seqs'
include { QIIME2_ALIGNMENT       } from '../modules/local/qiime2aln'
include { QIIME2_TREE            } from '../modules/local/qiime2tree'

// TODO: this can be captured in scripts and used in the versions.yml
include { SESSION_INFO           } from '../modules/local/rsessioninfo'

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

    def platform = params.platform.toLowerCase()

    if (!(["illumina","pacbio"].contains(platform))) {
        exit 1, "Only supported platforms (--platform argument) are currently 'pacbio' or 'illumina'"
    }

    // TODO: implement seqtable input
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

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // TODO: we may want to allow aggregation of the read files for larger projects;
    // current version is per read per sample
    PLOT_QUALITY_PROFILE (
        ch_samplesheet
    )

    ch_versions = ch_versions.mix(PLOT_QUALITY_PROFILE.out.versions.first())

    VSEARCH_EESTATS (
        ch_samplesheet
    )
    
    ch_versions = ch_versions.mix(VSEARCH_EESTATS.out.versions.first())

    // TODO: notice aggregation of data for multiqc and for version tracking, 
    // needs to be added throughout the workflow
    // ch_multiqc_files = ch_multiqc_files.mix(PLOTQUALITYPROFILE.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(PLOTQUALITYPROFILE.out.versions.first())

    // Subworkflows-Trimming and Filtering:
    //     cutadapt (overlapping paired: V4, COI)
    //     cutadapt (variable paired: ITS)
    //     cutadapt (long reads: PacBio 16S)
    //     DADA2 filterAndTrim
    //     Alternative est error filtering

    FILTER_AND_TRIM (
        ch_samplesheet
    )

    // Subworkflows-Denoising:
    //     DADA2 
    //          Pooled
    //          Per-sample
    //     Denoising (proposed alts): VSEARCH/Deblur
    //          Per-sample

    // TODO: Input for these should be the trimmed reads from above, but
    // they may need to be mixed in different ways depending on the 
    // denoising workflow used.  
    
    DADA2_DENOISE(FILTER_AND_TRIM.out.trimmed_infer)

    // Subworkflows-Taxonomic assignment (optional)
    ch_taxtab = Channel.empty()
    ch_boots =  Channel.empty()
    if (params.reference) {
        ref_file = file(params.reference, checkIfExists: true)
        species_file = params.species ? file(params.species, checkIfExists: true) : file("${projectDir}/assets/dummy_file")

        ASSIGN_TAXA_SPECIES(
            DADA2_DENOISE.out.readmap,
            ref_file,
            species_file
        )
        ch_taxtab = ASSIGN_TAXA_SPECIES.out.taxtab
        ch_boots = ASSIGN_TAXA_SPECIES.out.bootstraps
    }
    
    // Subworkflows-Alignment + Phylogenetic Tree (optional)
    DECIPHER(
        DADA2_DENOISE.out.nonchimeric_asvs
    )

    ch_tree = Channel.empty()

    // this seems like the sort of thing a function map 
    // would be useful for...
    // this needs to die with an error message if method is not defined, or
    // it needs to be caught above
    if (params.run_tree == 'phangorn') {
        PHANGORN(
            DECIPHER.out.alignment
        )
        ch_tree = PHANGORN.out.treeGTR
    } else if (params.run_tree == 'fasttree') {
        FASTTREE(
            DECIPHER.out.alignment
        )
        ch_tree = FASTTREE.out.treeGTR
    } 

    ROOT_TREE(
        ch_tree,
        params.run_tree
    )

    // QC
    READ_TRACKING(
        FILTER_AND_TRIM.out.trimmed_report,
        DADA2_DENOISE.out.seqtable_renamed,
        DADA2_DENOISE.out.inferred.collect(),
        DADA2_DENOISE.out.merged_seqs
    )

    PLOT_MERGED_HEATMAP(
        DADA2_DENOISE.out.merged_seqs
    )

    PLOT_ASV_DIST(
        DADA2_DENOISE.out.filtered_seqtable
    )

    // Subworkflow - Outputs

    SEQTABLE2TEXT(
        DADA2_DENOISE.out.filtered_seqtable
    )

    TAXTABLE2TEXT(
        ch_taxtab,
        ch_boots,
        DADA2_DENOISE.out.readmap
    )

    BIOM(
        DADA2_DENOISE.out.seqtable_renamed,
        ch_taxtab
    )

    QIIME2_FEATURETABLE(SEQTABLE2TEXT.out.seqtab2qiime)

    QIIME2_TAXTABLE(TAXTABLE2TEXT.out.taxtab2qiime)

    QIIME2_SEQUENCE(DADA2_DENOISE.out.nonchimeric_asvs)

    QIIME2_ALIGNMENT(DECIPHER.out.alignment)

    QIIME2_TREE(
        ch_tree,
        ROOT_TREE.out.rooted_tree
    )

    SESSION_INFO()

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
