#!/usr/bin/env nextflow

/*
========================================================================================
               D A D A 2   P I P E L I N E
========================================================================================
 DADA2 NEXTFLOW PIPELINE FOR UCT CBIO, HPCBio

----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    ===================================
     ${workflow.repository}/16S-rDNA-dada2-pipeline  ~  version ${params.version}
    ===================================
    Usage:

    This pipeline can be run specifying parameters in a config file or with command line flags.
    The typical example for running the pipeline with command line flags is as follows:
    nextflow run uct-cbio/16S-rDNA-dada2-pipeline --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex

    The typical command for running the pipeline with your own config (instead of command line flags) is as follows:
    nextflow run uct-cbio/16S-rDNA-dada2-pipeline -c dada2_user_input.config -profile uct_hex
    where:
    dada2_user_input.config is the configuration file (see example 'dada2_user_input.config')
    NB: -profile uct_hex still needs to be specified from the command line

    To override existing values from the command line, please type these parameters:

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. Currently profile available for UCT's HPC 'uct_hex' - create your own if necessary
                                    NB -profile should always be specified on the command line, not in the config file
      --trimFor                     integer. headcrop of read1 (set 0 if no trimming is needed)
      --trimRev                     integer. headcrop of read2 (set 0 if no trimming is needed)
      --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz)

    All available read preparation parameters:
      --trimFor                     integer. headcrop of read1
      --trimRev                     integer. headcrop of read2
      --truncFor                    integer. truncate read1 here (i.e. if you want to trim 10bp off the end of a 250bp R1, truncFor should be set to 240). enforced before trimFor/trimRev
      --truncRev                    integer. truncate read2 here ((i.e. if you want to trim 10bp off the end of a 250bp R2, truncRev should be set to 240). enforced before trimFor/trimRev
      --maxEEFor                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
      --maxEERev                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
      --truncQ                      integer. Truncate reads at the first instance of a quality score less than or equal to truncQ; default=2
      --maxN                        integer. Discard reads with more than maxN number of Ns in read; default=0
      --maxLen                      integer. maximum length of trimmed sequence; maxLen is enforced before trimming and truncation; default=Inf (no maximum)
      --minLen                      integer. minLen is enforced after trimming and truncation; default=50
      --rmPhiX                      {"T","F"}. remove PhiX from read
      --minOverlap                  integer. minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 integer. The maximum mismatches allowed in the overlap region; default=0
      --trimOverhang                {"T","F"}. If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off.
                                    "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction                                               primer region. Default="F" (false)

    Other arguments:
      --dadaOpt.XXX                 Set as e.g. --dadaOpt.HOMOPOLYMER_GAP_PENALTY=-1 Global defaults for the dada function, see ?setDadaOpt in R for available options and their defaults
      --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are
                                    pseudo pooling: "pseudo", true: "T", false: "F"
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run
                                    sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --idType                      The ASV IDs are renamed to simplify downstream analysis, in particular with downstream tools.  The
                                    default is "ASV", which simply renames the sequences in sequencial order.  Alternatively, this can be
                                    set to "md5" which will run MD5 on the sequence and generate a QIIME2-like unique hash.

    Help:
      --help                        Will print out summary above when executing nextflow run uct-cbio/16S-rDNA-dada2-pipeline

    Merging arguments (optional):
      --minOverlap                  The minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 The maximum mismatches allowed in the overlap region; default=0.
      --trimOverhang                If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off. "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen
                                    when reads are longer than the amplicon and read into the other-direction primer region. Default="F" (false)
      --minMergedLen                Minimum length of fragment *after* merging
      --maxMergedLen                Maximum length of fragment *after* merging

    Taxonomic arguments (optional):
      --species                     Specify path to fasta file. See dada2 addSpecies() for more detail.
    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// TODO: we need to validate/sanity-check more of the parameters
//Validate inputs

// ${deity} there has to be a better way to check this!
if ( (params.seqTables != false && params.reads != false) || 
     (params.input != false && params.reads != false) || 
     (params.seqTables != false && params.input != false )) {
    exit 1, "Only one of --reads, --input, or --seqTables is allowed!"
}

// Read-specific checks
if (params.reads != false) {
    if ( params.trimFor == false && params.amplicon == '16S') {
        exit 1, "Must set length of R1 (--trimFor) that needs to be trimmed (set 0 if no trimming is needed)"
    }

    if ( params.trimRev == false && params.amplicon == '16S') {
        exit 1, "Must set length of R2 (--trimRev) that needs to be trimmed (set 0 if no trimming is needed)"
    }

    // if ( params.reference == false ) {
    //     exit 1, "Must set reference database using --reference"
    // }

    if (params.fwdprimer == false && params.amplicon == 'ITS'){
        exit 1, "Must set forward primer using --fwdprimer"
    }

    if (params.revprimer == false && params.amplicon == 'ITS'){
        exit 1, "Must set reverse primer using --revprimer"
    }

    // TODO: add in single-end support here; will require an explicit flag
    Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { dada2ReadPairsToQual; dada2ReadPairsToDada2Qual; dada2ReadPairs }
} else if (params.seqTables != false) {
    // Experimental: combine pre-chimera sequence tables from prior DADA2 runs.  This assumes:
    // 1. Unique sample IDs in all runs
    // 2. Same DADA2 version used 
    // 3. Sequence table is the original DADA2 pre-chimera sequence table with the original sequence ASVs used as the ID

    Channel
        .fromFilePairs( params.seqTables, size: 1 )
        .ifEmpty { error "Cannot find any sequence tables matching: ${params.seqTables}" }
        .into { dada2SeqTabs }
} else if (params.input != false) {
    // Experimental: sample sheet input, per nf-core rules (CSV, with ID + FASTQ_R1 + FASTQ_R2)
    // This will get basic tuple-based support in place using code from the DSL2 rnaseq workflow:
    // https://github.com/nf-core/rnaseq/blob/e0dfce9af5c2299bcc2b8a74b6559ce055965455/subworkflows/local/input_check.nf
    // However this will become a separate subworkflow with DSL2 support later (see the link for DSL2 example)
    
    // see also: the check_samplesheet.py version in the bin directory, which does a high-level check on columns and data structure

    // process Check_SampleSheet {
    //     tag { "Check_SampleSheet" }
    //     // publishDir "${params.outdir}/FastQC-Raw", mode: "copy", overwrite: true

    //     input:
    //     file(samplesheet) from params.input

    //     output:
    //     file('final_samplesheet.csv') into samplesheet

    //     """
    //     check_samplesheet.py ${samplesheet} final_samplesheet.csv
    //     """        
    // }

    Channel
        .fromPath( params.input )
        .splitCsv(header:true, sep:',')
        .map{ row -> create_fastq_channel_simple(row) } // this doesn't really check files though...
        .into { dada2ReadPairsToQual; dada2ReadPairsToDada2Qual; dada2ReadPairs }
} else {
    exit 1, "Must set either --reads or --seqTables as input"
}

if (params.aligner == 'infernal' && params.infernalCM == false){
    exit 1, "Must set covariance model using --infernalCM when using Infernal"
}

if (!(['simple','md5'].contains(params.idType))) {
    exit 1, "--idType can currently only be set to 'simple' or 'md5', got ${params.idType}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

// Header log info
log.info "==================================="
log.info " ${params.base}/16S-rDNA-dada2-pipeline  ~  version ${params.version}"
log.info "==================================="
def summary = [:]
summary['Run Name']       = custom_runName ?: workflow.runName
summary['Reads']          = params.reads
summary['Forward primer'] = params.fwdprimer
summary['Reverse primer'] = params.revprimer
summary['Amplicon type']  = params.amplicon
summary['trimFor']        = params.trimFor
summary['trimRev']        = params.trimRev
summary['truncFor']       = params.truncFor
summary['truncRev']       = params.truncRev
summary['truncQ']         = params.truncQ
summary['maxEEFor']       = params.maxEEFor
summary['maxEERev']       = params.maxEERev
summary['maxN']           = params.maxN
summary['maxLen']         = params.maxLen
summary['minLen']         = params.minLen
summary['rmPhiX']         = params.rmPhiX
summary['minOverlap']     = params.minOverlap
summary['maxMismatch']    = params.maxMismatch
summary['trimOverhang']   = params.trimOverhang
summary['species']        = params.species
summary['dadaOpt']        = params.dadaOpt
summary['pool']           = params.pool
summary['qualityBinning'] = params.qualityBinning
summary['Reference']      = params.reference
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

/*
 *
 * Step 1: Filter and trim (run per sample?)
 *
 */

if (params.reads != false || params.input != false ) { // TODO maybe we should check the channel here

    process RunFastQC {
        tag { "FastQC-${pairId}" }
        publishDir "${params.outdir}/FastQC-Raw", mode: "copy", overwrite: true

        input:
        tuple val(pairId), file(in_fastq) from dada2ReadPairsToQual

        output:
        file '*_fastqc.{zip,html}' into fastqc_files

        """
        fastqc --nogroup -q ${in_fastq.get(0)} ${in_fastq.get(1)}
        """
    }

    process RunDADA2QC {
        tag { "DADA2-FASTQ-QC" }
        publishDir "${params.outdir}/dada2-RawQC", mode: "copy", overwrite: true

        input:
        path("fastq/*") from dada2ReadPairsToDada2Qual.flatMap({ n -> n[1] }).collect()

        when:
        params.skip_dadaQC == false

        output:
        file '*.pdf'

        script:
        template "DadaQC.R"
    }

    /* ITS amplicon filtering */

    // Note: should explore cutadapt options more: https://github.com/benjjneb/dada2/issues/785
    // https://cutadapt.readthedocs.io/en/stable/guide.html#more-than-one

    if (params.amplicon == 'ITS') {

        process ITSFilterAndTrimStep1 {
            tag { "ITS_Step1_${pairId}" }

            input:
            set pairId, file(reads) from dada2ReadPairs

            output:
            set val(pairId), "${pairId}.R[12].noN.fastq.gz" optional true into itsStep2
            set val(pairId), "${pairId}.out.RDS" into itsStep3Trimming  // needed for join() later
            file('forward_rc') into forwardP
            file('reverse_rc') into reverseP

            when:
            params.precheck == false

            script:
            template "ITSFilterAndTrimStep1.R"
        }
        
        process ITSFilterAndTrimStep2 {
            tag { "ITS_Step2_${pairId}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            set pairId, reads from itsStep2
            file(forP) from forwardP
            file(revP) from reverseP
            
            output:
            set val(pairId), "${pairId}.R[12].cutadapt.fastq.gz" optional true into itsStep3
            file "*.cutadapt.out" into cutadaptToMultiQC

            when:
            params.precheck == false

            script:
            """
            FWD_PRIMER=\$(<forward_rc)
            REV_PRIMER=\$(<reverse_rc)
            
            cutadapt -g "${params.fwdprimer}" -a \$FWD_PRIMER \\
                -G "${params.revprimer}" -A \$REV_PRIMER \\
                --cores ${task.cpus} \\
                -n 2 \\
                -o "${pairId}.R1.cutadapt.fastq.gz" \\
                -p "${pairId}.R2.cutadapt.fastq.gz" \\
                "${reads[0]}" "${reads[1]}" > "${pairId}.cutadapt.out"
            """
        }

        process ITSFilterAndTrimStep3 {
            tag { "ITS_Step3_${pairId}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            set pairId, file(reads), file(trimming) from itsStep3.join(itsStep3Trimming)

            output:
            set val(pairId), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" optional true into filteredReadsforQC, filteredReads
            tuple val("R1"), file("*.R1.filtered.fastq.gz") optional true into forReadsLE
            tuple val("R2"), file("*.R2.filtered.fastq.gz") optional true into revReadsLE
            file "*.trimmed.txt" into trimTracking

            when:
            params.precheck == false

            script:
            template "ITSFilterAndTrimStep3.R"
        }
        
    }
    /* 16S amplicon filtering */
    else if (params.amplicon == '16S'){
        process FilterAndTrim {
            tag { "FilterAndTrim_${pairId}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            tuple pairId, file(reads) from dada2ReadPairs

            output:
            tuple val(pairId), "${pairId}.R1.filtered.fastq.gz", "${pairId}.R2.filtered.fastq.gz" optional true into filteredReadsforQC,filteredReads
            tuple val("R1"), file("${pairId}.R1.filtered.fastq.gz") optional true into forReadsLE
            tuple val("R2"), file("${pairId}.R2.filtered.fastq.gz") optional true into revReadsLE
            file "*.trimmed.txt" into trimTracking

            when:
            params.precheck == false

            script:
            phix = params.rmPhiX ? '--rmPhiX TRUE' : '--rmPhiX FALSE'
            template "FilterAndTrim.R"
        }
        cutadaptToMultiQC = Channel.empty()
    } else {
        // We need to shut this down!
        Channel.empty().into {cutadaptToMultiQC;filteredReads;filteredReadsforQC}
    }

    process RunFastQC_postfilterandtrim {
        tag { "FastQC_post_FT.${pairId}" }
        publishDir "${params.outdir}/FastQC-Post-FilterTrim", mode: "copy", overwrite: true

        input:
        set val(pairId), file(filtFor), file(filtRev) from filteredReadsforQC

        output:
        file '*_fastqc.{zip,html}' into fastqc_files_post

        when:
        params.precheck == false

        """
        fastqc --nogroup -q ${filtFor} ${filtRev}
        """
    }

    process RunMultiQC_postfilterandtrim {
        tag { "runMultiQC_postfilterandtrim" }
        publishDir "${params.outdir}/MultiQC-Post-FilterTrim", mode: 'copy', overwrite: true

        input:
        file('./RawSeq/*') from fastqc_files.collect().ifEmpty([])
        file('./TrimmedSeq/*') from fastqc_files_post.collect().ifEmpty([])
        file('./Cutadapt/*') from cutadaptToMultiQC.collect().ifEmpty([])

        output:
        file "*_report.html" into multiqc_report_post
        file "*_data"

        when:
        params.precheck == false & params.skip_multiQC == false

        script:
        interactivePlots = params.interactiveMultiQC == true ? "-ip" : ""
        """
        multiqc ${interactivePlots} .
        """
    }

    process MergeTrimmedTable {
        tag { "MergeTrimmedTable" }
        publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

        input:
        file trimData from trimTracking.collect()

        output:
        file "all.trimmed.csv" into trimmedReadTracking

        when:
        params.precheck == false

        script:
        template "MergeTrimTable.R"
    }

    /*
     *
     * Step 2: Learn error rates (run on all samples)
     * 
     * Note: this step has been re-configured to run R1 and R2 batches in parallel
     *
     */

    // the incoming data needs to be processed per R1/R2; this next step assists with this
    forReadsLE
            .concat(revReadsLE)
            .groupTuple(sort: true)
            .into { ReadsLE;ReadsInfer;ReadsMerge }

    process LearnErrors {
        tag { "LearnErrors:${readmode}" }
        publishDir "${params.outdir}/dada2-LearnErrors", mode: "copy", overwrite: true

        input:
        tuple val(readmode), file(reads) from ReadsLE

        output:
        tuple val(readmode), file("errors.${readmode}.RDS") into errorModels
        file "${readmode}.err.pdf"

        when:
        params.precheck == false

        script:
        dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
        template "LearnErrors.R"
    }

    /*
     *
     * Step 3: Sample Inference, Merge Pairs
     *
     */

    /*
     *
     * Step 4: Construct sequence table
     *
     */

    if (params.pool == "T" || params.pool == 'pseudo') {

        process DadaInfer {
            tag { "DadaInfer:${readmode}" }
            publishDir "${params.outdir}/dada2-Derep-Pooled", mode: "copy", overwrite: true

            input:
            tuple val(readmode), file(err), file(reads) from errorModels
                .join(ReadsInfer)

            output:
            // Note that the mode ('merged', 'R1', 'R2') can now potentially allow SE read analysis
            tuple val(readmode), file("all.dd.${readmode}.RDS") into dadaMerge,dadaToReadTracking,dadaSECalls

            when:
            params.precheck == false

            script:
            dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
            template "DadaPooled.R"
        }

        // Merge the pooled runs
        process MergePooled {
            tag { "mergeDadaRDS" }
            publishDir "${params.outdir}/dada2-Derep-Pooled", mode: "copy", overwrite: true

            input:
            file(dds) from dadaMerge
                .map { it[1] }
                .collect()
            file(filts) from ReadsMerge
                .map { it[1] }
                .flatten()
                .collect()

            output:
            tuple val("merged"), file("seqtab.merged.RDS") into seqTable,rawSeqTableToRename
            file "all.mergers.RDS" into mergerTracking
            file "seqtab.original.merged.RDS" // we keep this for comparison and possible QC

            when:
            params.precheck == false

            script:
            template "MergePairs.R"
        }

    } else {
        // pool = F, process per sample or in batches (a little more compute efficient)

        process PerSampleInferDerepAndMerge {
            tag { "PerSampleInferDerepAndMerge" }
            publishDir "${params.outdir}/dada2-Derep-Single/Per-Sample", mode: "copy", overwrite: true

            input:
            tuple val(pairId), file(r1), file(r2) from filteredReads
            file errs from errorModels.collect()

            output:
            tuple val(pairId), file("${pairId}.merged.RDS") into mergedReads
            tuple val(pairId), file("${pairId}.dd.R{1,2}.RDS") into perSampleDadaToMerge

            when:
            params.precheck == false

            script:
            dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
            template "PerSampleDadaInfer.R"
        }

        process MergeDadaRDS {
            tag { "mergeDadaRDS" }
            publishDir "${params.outdir}/dada2-Derep-Single", mode: "copy", overwrite: true

            input:
            file(dds)from perSampleDadaToMerge
                          .map { it[1] }
                          .flatten()
                          .collect()

            output:
            tuple val(readmode), file("all.dd.R{1,2}.RDS") into dadaToReadTracking,dadaSECalls

            when:
            params.precheck == false

            script:
            template "MergePerSampleDada.R"
        }

        process SequenceTable {
            tag { "SequenceTable" }
            publishDir "${params.outdir}/dada2-Derep-Single", mode: "copy", overwrite: true

            input:
            file mr from mergedReads.collect()

            output:
            tuple val("merged"), file("seqtab.merged.RDS") into seqTable,rawSeqTableToRename
            file "all.mergers.RDS" into mergerTracking
            file "seqtab.original.merged.RDS" // we keep this for comparison and possible QC
            
            when:
            params.precheck == false

            script:
            template "PerSampleSeqTable.R"
        }
    }

    if (params.processSE) {
        process SESequenceTable {
            tag { "SESequenceTable:${readmode}" }
            publishDir "${params.outdir}/dada2-SE", mode: "copy", overwrite: true

            input:
            tuple val(readmode), file(dds) from dadaSECalls

            output:
            tuple val(readmode), file("seqtab.${readmode}.RDS") into SEChimera,RawSEChimeraToRename

            script:
            """
            #!/usr/bin/env Rscript
            suppressPackageStartupMessages(library(dada2))
            dd <- readRDS("${dds}")
            seqtab <- makeSequenceTable(dd)
            saveRDS(seqtab, "seqtab.${readmode}.RDS")
            """
        }

        // add to the queue with the merged sequences

    } else {
        Channel.empty().into { SEChimera;RawSEChimeraToRename }
    }
} else if (params.seqTables != false) { // TODO maybe we should check the channel here
    process MergeSeqTables {
        tag { "MergeSeqTables" }
        publishDir "${params.outdir}/dada2-MergedSeqTable", mode: 'copy'

        input:
        file(st) from dada2SeqTabs
                    .map { it[1] }
                    .collect()

        output:
        tuple val("merged"), file("seqtab.merged.RDS") into seqTable, rawSeqTableToRename

        script:
        template "MergeSeqTables.R"
    }
    Channel.empty().into { SEChimera;RawSEChimeraToRename;trimmedReadTracking;dadaToReadTracking;mergerTracking }
} 

/*
 *
 * Step 8: Remove chimeras
 *
 */

if (!params.skipChimeraDetection) {
    process RemoveChimeras {
        tag { "RemoveChimeras:${seqtype}" }
        publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-Chimera" : "${params.outdir}/dada2-Chimera/${seqtype}"
            }, mode: 'copy'

        input:
        tuple val(seqtype), file(st) from seqTable.concat(SEChimera)

        output:
        tuple val(seqtype), file("seqtab_final.${seqtype}.RDS") into seqTableToTax,seqTableToRename

        when:
        params.precheck == false

        script:
        chimOpts = params.removeBimeraDenovoOptions != false ? ", ${params.removeBimeraDenovoOptions}" : ''
        template "RemoveChimeras.R"
    }
} else {
    if (params.skipChimeraDetection == true) {
        seqTable.into {seqTableToTax;seqTableToRename}
    } else if (params.skipChimeraDetection == 'all') { // stop completely!
        // this should effectively halt the pipeline from going further
        Channel.empty()
            .into {seqTableToTax;seqTableToRename}
    }
}

/*
 *
 * Step 9: Taxonomic assignment
 *
 */

if (params.reference) {

    if (params.taxassignment == 'rdp') {
        // TODO: we could combine these into the same script
        refFile = file(params.reference)

        if (params.species) {

            speciesFile = file(params.species)

            process AssignTaxSpeciesRDP {
                tag { "AssignTaxSpeciesRDP:${seqtype}" }
                publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Taxonomy" : "${params.outdir}/dada2-Taxonomy/${seqtype}"
                    }, mode: "copy"

                input:
                tuple val(seqtype), file(st) from seqTableToTax
                file ref from refFile
                file sp from speciesFile

                output:
                tuple val(seqtype), file("tax_final.${seqtype}.RDS") into taxFinal,taxTableToTable
                tuple val(seqtype), file("bootstrap_final.${seqtype}.RDS") into bootstrapFinal

                when:
                params.precheck == false

                script:
                template "AssignTaxonomySpeciesRDP.R"
            }

        } else {

            process AssignTaxonomyRDP {
                tag { "TaxonomyRDP:${seqtype}" }
                publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Taxonomy" : "${params.outdir}/dada2-Taxonomy/${seqtype}"
                    }, mode: "copy", overwrite: true

                input:
                tuple val(seqtype), file(st) from seqTableToTax
                file ref from refFile

                output:
                tuple val(seqtype), file("tax_final.${seqtype}.RDS") into taxFinal,taxTableToTable
                tuple val(seqtype), file("bootstrap_final.${seqtype}.RDS") into bootstrapFinal

                when:
                params.precheck == false

                // TODO: this is not tested yet
                script:
                taxLevels = params.taxLevels ? "c( ${params.taxLevels} )," : ''
                template "AssignTaxonomyRDP.R"
            }
        }
    } else if (params.taxassignment == 'idtaxa') {
        // Experimental!!! This assigns full taxonomy to species level, but only for
        // some databases; unknown whether this works with concat sequences.  ITS
        // doesn't seem to be currently supported
        process TaxonomyIDTAXA {
            tag { "TaxonomyIDTAXA:${seqtype}" }
            publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Taxonomy" : "${params.outdir}/dada2-Taxonomy/${seqtype}"
                    }, mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(st) from seqTableToTax
            file ref from refFile // this needs to be a database from the IDTAXA site

            output:
            tuple val(seqtype), file("tax_final.${seqtype}.RDS") into taxFinal,taxTableToTable
            tuple val(seqtype), file("bootstrap_final.${seqtype}.RDS") into bootstrapFinal
            file "raw_idtaxa.${seqtype}.RDS"

            when:
            params.precheck == false

            // TODO: this is not tested yet
            script:
            template "TaxonomyIDTAXA.R"
        }

    } else if (params.taxassignment) {
        exit 1, "Unknown taxonomic assignment method set: ${params.taxassignment}"
    } else {
        exit 1, "No taxonomic assignment method set, but reference passed"
    }
} else {
    // set tax channels to 'false', do NOT assign taxonomy
    Channel.empty().into { taxFinal;taxTableToTable;bootstrapFinal }
}

// Note: this is currently a text dump.  We've found the primary issue with
// downstream analysis is getting the data in a form that can be useful as
// input, and there isn't much consistency with this as of yet.  So for now
// we're using the spaghetti approach (see what sticks).  Also, we  are running
// into issues with longer sequences (e.g. concatenated ones) used as IDs with
// tools like Fasttree (it doesn't seem to like that).

// Safest way may be to save the simpleID -> seqs as a mapping file, use that in
// any downstream steps (e.g. alignment/tree), then munge the seq names back
// from the mapping table

/*
 *
 * Step 8.5: Rename ASVs
 *
 * A number of downstream programs have issues with sequences as IDs, here we
 * (optionally) rename these
 *
 */

process RenameASVs {
    tag { "RenameASVs:${seqtype}" }
    publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Tables" : "${params.outdir}/dada2-Tables/${seqtype}"
                    }, mode: "copy", overwrite: true

    input:
    // This is a tricky join; we have two incoming channels that we want to
    // join by sequence type (merged, R1, R2). seqTableToRename can have all
    // three, but rawSeqTableToRename (originating from merging) has only
    // one, while SE-specific seqtables (with RawSEChimeraToRename) are 
    // generated on demand later, primarily to act as a switch.  So we join() the 
    // channels by ID, but we concatenate the original seqtable channels together 

    tuple val(seqtype), file(st), file(rawst) from seqTableToRename
        .join(rawSeqTableToRename
                .concat(RawSEChimeraToRename))

    output:
    tuple val(seqtype), file("seqtab_final.${params.idType}.${seqtype}.RDS") into seqTableFinalToBiom,seqTableFinalToTax,seqTableFinalTree,seqTableFinalTracking,seqTableToTable,seqtabToPhyloseq,seqtabToTaxTable
    tuple val(seqtype), file("asvs.${params.idType}.nochim.${seqtype}.fna") into seqsToAln, seqsToQIIME2
    tuple val(seqtype), file("readmap.${seqtype}.RDS") into readsToRenameTaxIDs // needed for remapping tax IDs
    file "asvs.${params.idType}.raw.${seqtype}.fna"

    script:
    template "RenameASVs.R"
}

process GenerateSeqTables {
    tag { "GenerateSeqTables:${seqtype}" }
    publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Tables" : "${params.outdir}/dada2-Tables/${seqtype}"
                    }, mode: "copy", overwrite: true

    input:
    tuple val(seqtype), file(st) from seqTableToTable

    output:
    tuple val(seqtype), file("seqtab_final.${params.idType}.${seqtype}.qiime2.txt") into featuretableToQIIME2
    file "*.txt"

    when:
    params.precheck == false

    script:
    template "GenerateSeqTables.R"
}

process GenerateTaxTables {
    tag { "GenerateTaxTables:${seqtype}" }
    publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Tables" : "${params.outdir}/dada2-Tables/${seqtype}"
                    }, mode: "copy", overwrite: true

    input:
    tuple val(seqtype), file(tax), file(bt), file(map) from taxTableToTable.join(bootstrapFinal).join(readsToRenameTaxIDs)

    output:
    tuple val(seqtype), file("tax_final.${params.idType}.${seqtype}.RDS") into taxtabToPhyloseq
    tuple val(seqtype), file("tax_final.${params.idType}.${seqtype}.txt") into taxtableToQIIME2
    file "*.txt"

    when:
    params.precheck == false

    script:
    template "GenerateTaxTables.R"
}

/*
 *
 * Step 10: Align and construct phylogenetic tree
 *
 */


/*
 *
 * Step 10a: Alignment
 *
 */

// NOTE: 'when' directive doesn't work if channels have the same name in
// two processes

if (!params.precheck && params.runTree && params.amplicon != 'ITS') {

    if (params.aligner == 'infernal') {

        cmFile = file(params.infernalCM)

        process AlignReadsInfernal {
            tag { "AlignReadsInfernal:${seqtype}" }
            publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Infernal" : "${params.outdir}/dada2-Infernal/${seqtype}"
                    }, mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(seqs) from seqsToAln
            file cm from cmFile

            output:
            file "aligned_seqs.${seqtype}.stk"
            file "aln.${seqtype}.scores"
            tuple val(seqtype), file("aligned_seqs.${seqtype}.fasta") optional true into alnFile,alnToQIIME2

            script:
            """
            # from the original IM-TORNADO pipeline
            cmalign --cpu ${task.cpus} \\
                  -g --notrunc --sub --dnaout --noprob \\
                  --sfile aln.${seqtype}.scores \\
                  -o aligned_seqs.${seqtype}.stk \\
                  ${cm} ${seqs}

            # script from P. Jeraldo (Mayo)
            stkToFasta.py aligned_seqs.${seqtype}.stk aligned_seqs.${seqtype}.fasta
            """
        }
    } else if (params.aligner == 'DECIPHER') {

        process AlignReadsDECIPHER {
            tag { "AlignReadsDECIPHER:${seqtype}" }
            publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-DECIPHER" : "${params.outdir}/dada2-DECIPHER/${seqtype}"
                    }, mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(seqs) from seqsToAln

            output:
            tuple val(seqtype), file("aligned_seqs.${seqtype}.fasta") optional true into alnFile,alnToQIIME2
            
            script:
            template "AlignReadsDECIPHER.R"
        }
    } else {
        exit 1, "Unknown aligner option: ${params.aligner}"
    }

    /*
     *
     * Step 10b: Construct phylogenetic tree
     *
     */

    if (params.runTree == 'phangorn') {

        process GenerateTreePhangorn {
            tag { "GenerateTreePhangorn:${seqtype}" }
            publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Phangorn" : "${params.outdir}/dada2-Phangorn/${seqtype}"
                    }, mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(aln) from alnFile

            output:
            file "phangorn.tree.${seqtype}.RDS" into treeRDS
            file "tree.${seqtype}.newick" into treeFile
            file "tree.GTR.${seqtype}.newick" into treeGTRFile

            script:
            template "PhangornML.R"
        }
    } else if (params.runTree == 'fasttree') {

        process GenerateTreeFasttree {
            tag { "GenerateTreeFasttree:${seqtype}" }
            publishDir path: {
                    seqtype == "merged" ? "${params.outdir}/dada2-Fasttree" : "${params.outdir}/dada2-Fasttree/${seqtype}"
                    }, mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(aln) from alnFile

            output:
            tuple val(seqtype), file("fasttree.${seqtype}.tree") into treeGTRFile, treeToQIIME2
            // need to deadend the other channels, they're hanging here

            script:
            """
            OMP_NUM_THREADS=${task.cpus} FastTree -nt \\
                -gtr -gamma -spr 4 -mlacc 2 -slownni \\
                -out fasttree.${seqtype}.tree \\
                ${aln}
            """
        }

    } else {
        // dead-end channels generated above
    }

    process RootTree {
        tag { "RootTree:${seqtype}" }
        publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-RootedTree" : "${params.outdir}/dada2-RootedTree/${seqtype}"
            }, mode: "copy"

        input:
        tuple val(seqtype), file(tree) from treeGTRFile

        output:
        tuple val(seqtype), file("rooted.${seqtype}.newick") into rootedTreeFile, rootedToQIIME2
        // need to deadend the other channels, they're hanging here

        script:
        template "RootTree.R"
    }
} else {
    // Note these are caught downstream
    Channel.empty().into { alnToQIIME2;treeToQIIME2;rootedToQIIME2 }
}

// TODO: rewrite using the python BIOM tools?

process ToBiomFile {
    tag { "ToBiomFile:${seqtype}" }
    publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-BIOM" : "${params.outdir}/dada2-BIOM/${seqtype}"
            }, mode: "copy", overwrite: true

    input:
    tuple val(seqtype), file(sTable), file(tTable) from seqTableFinalToBiom.join(taxFinal)

    output:
    file "dada2.${seqtype}.biom" into biomFile

    when:
    params.toBIOM == true

    script:
    template "ToBIOM.R"
}

/*
 *
 * Step 10: Track reads
 *
 */

process ReadTracking {
    tag { "ReadTracking" }
    publishDir "${params.outdir}/dada2-ReadTracking", mode: "copy", overwrite: true

    input:
    file trimmedTable from trimmedReadTracking
    tuple(seqtype), file(sTable) from seqTableFinalTracking
    file dds from dadaToReadTracking.collect()
    file mergers from mergerTracking

    output:
    file "all.readtracking.${seqtype}.txt"

    script:
    template "ReadTracking.R"
}

if (params.toQIIME2) {

    process ToQIIME2FeatureTable {
        tag { "QIIME2-SeqTable:${seqtype}" }
        label 'QIIME2'
        publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-QIIME2" : "${params.outdir}/dada2-QIIME2/${seqtype}"
            }, mode: "copy"

        input:
        tuple val(seqtype), file(seqtab) from featuretableToQIIME2

        output:
        file "*.qza"

        script:
        """
        biom convert -i ${seqtab} \
            -o seqtab-biom-table.biom \
            --table-type="OTU table" \
            --to-hdf5

        qiime tools import \
            --input-path seqtab-biom-table.biom \
            --input-format BIOMV210Format \
            --output-path seqtab_final.${params.idType}.${seqtype}.qza \
            --type 'FeatureTable[Frequency]'
        """
    }

    process ToQIIME2TaxTable {
        tag { "QIIME2-Taxtable:${seqtype}" }
        label 'QIIME2'
        publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-QIIME2" : "${params.outdir}/dada2-QIIME2/${seqtype}"
            }, mode: "copy"

        input:
        tuple val(seqtype), file(taxtab) from taxtableToQIIME2

        output:
        file "*.qza"

        when:
        taxtableToQIIME2 != false

        script:
        """
        tail -n +2 ${taxtab} > headerless.txt
        qiime tools import \
            --input-path headerless.txt \
            --input-format HeaderlessTSVTaxonomyFormat \
            --output-path tax_final.${params.idType}.${seqtype}.qza \
            --type 'FeatureData[Taxonomy]'
        """
    }

    process ToQIIME2Seq {
        tag { "QIIME2-Seq:${seqtype}" }
        label 'QIIME2'
        publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-QIIME2" : "${params.outdir}/dada2-QIIME2/${seqtype}"
            }, mode: "copy"

        input:
        tuple val(seqtype), file(seqs) from seqsToQIIME2

        output:
        file "*.qza"

        script:
        """
        qiime tools import \
            --input-path ${seqs} \
            --output-path sequences.${seqtype}.qza \
            --type 'FeatureData[Sequence]'
        """
    }

    process ToQIIME2Aln {
        tag { "QIIME2-Aln:${seqtype}" }
        label 'QIIME2'
        publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-QIIME2" : "${params.outdir}/dada2-QIIME2/${seqtype}"
            }, mode: "copy"

        input:
        tuple val(seqtype), file(aln) from alnToQIIME2

        when:
        alnToQIIME2 != false

        output:
        file "*.qza"

        script:
        """
        qiime tools import \
            --input-path ${aln} \
            --output-path aligned-sequences.${seqtype}.qza \
            --type 'FeatureData[AlignedSequence]'
        """
    }

    process ToQIIME2Tree {
        tag { "QIIME2-Tree:${seqtype}" }
        label 'QIIME2'
        publishDir path: {
            seqtype == "merged" ? "${params.outdir}/dada2-QIIME2" : "${params.outdir}/dada2-QIIME2/${seqtype}"
            }, mode: "copy"

        input:
        tuple val(seqtype), file(rooted), file(tree) from rootedToQIIME2.join(treeToQIIME2)

        output:
        file "*.qza"

        script:
        """
        qiime tools import \
            --input-path ${tree} \
            --output-path unrooted-tree.${seqtype}qza \
            --type 'Phylogeny[Unrooted]'

        qiime tools import \
            --input-path ${rooted} \
            --output-path rooted-tree.${seqtype}.qza \
            --type 'Phylogeny[Rooted]'
        """
    }
}

//TODO: this could eventually go into a report process to consolidate workflow info;
//      deciding between pander, knitr, or combination of the two

process SessionInfo {
        tag { "R-sessionInfo" }
        label 'sessionInfo'
        publishDir "${params.outdir}/SessionInfo", mode: "link"

        output:
        file "sessionInfo.Rmd"

        script:
        template "SessionInfo.R"
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    def subject = "[${params.base}/16S-rDNA-dada2-pipeline] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[${params.base}/16S-rDNA-dada2-pipeline] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    // if (params.email) {
    //     try {
    //       if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
    //       // Try to send HTML e-mail using sendmail
    //       [ 'sendmail', '-t' ].execute() << sendmail_html
    //       log.info "[${params.base}/16S-rDNA-dada2-pipeline] Sent summary e-mail to $params.email (sendmail)"
    //     } catch (all) {
    //       // Catch failures and try with plaintext
    //       [ 'mail', '-s', subject, params.email ].execute() << email_txt
    //       log.info "[${params.base}/16S-rDNA-dada2-pipeline] Sent summary e-mail to $params.email (mail)"
    //     }
    // }
}


// code modified from the nf-core RNA-Seq workflow
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def sample_data = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        sample_data = [ row.sample, file(row.fastq_1) ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        sample_data = [ row.sample, file(row.fastq_1), file(row.fastq_2) ] 
    }
    return sample_data
}

def create_fastq_channel_simple(LinkedHashMap row) {
    // add path(s) of the fastq file(s) to the meta map
    def sample_data = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    sample_data = [ row.sample, [ file(row.fastq_1), file(row.fastq_2) ] ] 
    return sample_data
}
