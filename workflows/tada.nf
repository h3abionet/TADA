/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { PLOTQUALITYPROFILE     } from '../modules/local/plotqualityprofile'
// TODO: we will use the module, but the trimming step will get more complex and
//       should be moved to a subworkflow
// include { FILTERANDTRIM          } from '../subworkflows/local/'
include { FILTERANDTRIM          } from '../modules/local/filterandtrim'
include { MERGETRIMTABLES        } from '../modules/local/mergetrimtables'
include { LEARNERRORS            } from '../modules/local/learnerrors'
include { DADAINFER              } from '../modules/local/dadainfer'
include { POOLEDSEQTABLE         } from '../modules/local/pooledseqtable'
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

    // def platform = ''
    // platform = params.platform.toLowerCase()

    // if (!(["illumina","pacbio","pacbio-kinnex"].contains(platform))) {
    //     exit 1, "Only supported platforms (--platform argument) are currently 'pacbio', 'pacbio-kinnex', or 'illumina'"
    // }

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

    if (params.aligner == 'infernal' && params.infernalCM == false){
        exit 1, "Must set covariance model using --infernalCM when using Infernal"
    }

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
    PLOTQUALITYPROFILE (
        ch_samplesheet
    )
    
    // ch_multiqc_files = ch_multiqc_files.mix(PLOTQUALITYPROFILE.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(PLOTQUALITYPROFILE.out.versions.first())

    // Subworkflows-Trimming and Filtering:
    //     cutadapt (overlapping paired: V4, COI)
    //     cutadapt (variable paired: ITS)
    //     cutadapt (long reads: PacBio 16S)
    //     DADA2 filterAndTrim
    //     Alternative est error filtering

    FILTERANDTRIM (
        ch_samplesheet
    )

    ch_trimmed_reads = FILTERANDTRIM.out.trimmed
    ch_reports = FILTERANDTRIM.out.trimmed_report.collect()

    
    MERGETRIMTABLES(
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
    ch_trimmed_infer = FILTERANDTRIM.out.trimmed_R1
            // .concat(filteredReadsR2.ifEmpty([]))
            .map { [ 'R1', it[1]] }
            .concat(FILTERANDTRIM.out.trimmed_R2.map {['R2', it[1]] } )
            .groupTuple(sort: true)

    // Subworkflows-Denoising:
    //     DADA2 
    //          Pooled
    //          Per-sample
    //     Denoising (proposed alt): VSEARCH/Deblur
    //          Per-sample

    // TODO: This can go into a DADA2-specific subworkflow
    LEARNERRORS (
        ch_trimmed_infer
    )

    ch_infer = LEARNERRORS.out.error_models.join(ch_trimmed_infer)

    // this is always in pooled mode at the moment, should be adjusted
    DADAINFER(
        ch_infer
    )

    ch_trimmed = ch_trimmed_infer
        .map { it[1] }
        .flatten()
        .collect()

    POOLEDSEQTABLE(
        DADAINFER.out.inferred.collect(),
        ch_trimmed
        )
    // if (params.pool == "T" || params.pool == 'pseudo') { 

        // process DadaInfer {
        //     tag { "DadaInfer:${readmode}" }
        //     publishDir "${params.outdir}/dada2-Derep-Pooled", mode: "copy", overwrite: true

        //     input:
        //     // DADA2 step runs on all R1 and/or on all R2
        //     tuple val(readmode), file(err), file(reads) from errorModelsPooled
        //         .join(ReadsInfer)

        //     output:
        //     // Note that the mode ('merged', 'R1', 'R2') can now potentially allow SE read analysis
        //     file("all.dd.${readmode}.RDS") into dadaMerge,dadaToReadTracking

        //     when:
        //     params.precheck == false

        //     script:
        //     dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
        //     template "DadaPooled.R"
        // }

        // // This one is a little tricky. We can't know a priori how many instances of reads (R1 and R2) 
        // // are present outside the process, but we can determine this internally within the process 
        // // when we collect all of them.
        // // So, here we check the size of the collected channel containing the denoised models; if 
        // // there are two then this is a paired-end run, otherwise it's single-end. Logic is in the R script
        
        // process PooledSeqTable {
        //     tag { "PooledSeqTable:${readmode}" }
        //     publishDir "${params.outdir}/dada2-OriginalSeqTable", mode: "copy", overwrite: true

        //     input:
        //     // we don't care about the mode here, so only get the dds (dada-inferred) RDS files
        //     file(dds) from dadaMerge.collect()
        //     // we don't care about the mode here, we only grab the reads
        //     file(filts) from ReadsMerge
        //         .map { it[1] }
        //         .flatten()
        //         .collect()

        //     output:
        //     tuple val(readmode), file("seqtab.${readmode}.RDS") into seqTable,rawSeqTableToRename
        //     file "all.merged.RDS" optional true into mergerTracking,mergerQC
        //     file "seqtab.original.*.RDS" into seqtabQC// we keep this for comparison and possible QC

        //     when:
        //     params.precheck == false

        //     script:
        //     // We could switch this to 'paired' vs 'single-end' as well
        //     readmode = dds.size() == 2 ? 'merged' : 'R1'
        //     template "SeqTables.R"
        // }
    // } else {

        // process PerSampleInferDerepAndMerge {
        //     tag { "PerSampleInferDerepAndMerge:${meta.id}" }
        //     publishDir "${params.outdir}/dada2-Derep-Single/Per-Sample", mode: "copy", overwrite: true

        //     input:
        //     tuple val(meta), file(reads) from readsToPerSample
        //     file(errs) from errorModelsPerSample.collect()

        //     output:
        //     file("${meta.id}.{R1,merged}.RDS") into combinedReads
        //     tuple val(meta), file("${meta.id}.dd.R{1,2}.RDS") into perSampleDadaToMerge
        //     val(readmode) into modeSeqTable

        //     when:
        //     params.precheck == false

        //     script:
        //     dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
        //     readmode = errs.size() == 2 ? 'merged' : 'R1'
        //     template "PerSampleDadaInfer.R"
        // }

        // process MergeDadaRDS {
        //     tag { "mergeDadaRDS" }
        //     publishDir "${params.outdir}/dada2-Derep-Single", mode: "copy", overwrite: true

        //     input:
        //     file(dds) from perSampleDadaToMerge
        //                   .map { it[1] }
        //                   .flatten()
        //                   .collect()

        //     output:
        //     file("all.dd.R{1,2}.RDS") into dadaToReadTracking

        //     when:
        //     params.precheck == false

        //     script:
        //     template "MergePerSampleDada.R"
        // }

        // process SequenceTable {
        //     tag { "SequenceTable:${readmode}" }
        //     publishDir "${params.outdir}/dada2-Derep-Single", mode: "copy", overwrite: true

        //     input:
        //     file(mr) from combinedReads.collect()
        //     val(readmode) from modeSeqTable.first()

        //     output:
        //     tuple val(readmode), file("seqtab.${readmode}.RDS") into seqTable,rawSeqTableToRename
        //     file "all.merged.RDS" optional true into mergerTracking,mergerQC
        //     file "seqtab.original.${readmode}.RDS" into seqtabQC // we keep this for comparison and possible QC
            
        //     when:
        //     params.precheck == false

        //     script:
        //     template "PerSampleSeqTable.R"
        // }
    // }

    // } else if (params.seqTables) { // TODO maybe we should check the channel here
    //     process MergeSeqTables {
    //         tag { "MergeSeqTables" }
    //         publishDir "${params.outdir}/dada2-MergedSeqTable", mode: 'copy'

    //         input:
    //         file(st) from dada2SeqTabs
    //                     .map { it[1] }
    //                     .collect()

    //         output:
    //         tuple val("merged"), file("seqtab.merged.RDS") into seqTable, rawSeqTableToRename

    //         script:
    //         template "MergeSeqTables.R"
    //     }
    //     Channel.empty().into { SEChimera;RawSEChimeraToRename;trimmedReadTracking;dadaToReadTracking;mergerTracking;mergerQC }
    // }    

    // Subworkflows-Taxonomic assignment (optional)

    // Subworkflows-Alignment + Phylogenetic Tree

    // Subworkflows-Alternative outputs

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
