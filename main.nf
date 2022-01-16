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

// TODO: add checks on options

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

//Validate inputs
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

if (params.aligner == 'infernal' && params.infernalCM == false){
    exit 1, "Must set covariance model using --infernalCM when using Infernal"
}

if (!(['simple','md5'].contains(params.idType))) {
    exit 1, "--idType can only be set to 'simple' or 'md5', got ${params.idType}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { dada2ReadPairsToQual; dada2ReadPairsToDada2Qual; dada2ReadPairs }

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

process runFastQC {
    tag { "FASTQC-${pairId}" }
    publishDir "${params.outdir}/FASTQC-Raw", mode: "copy", overwrite: true

    input:
    set pairId, file(in_fastq) from dada2ReadPairsToQual

    output:
    file '*_fastqc.{zip,html}' into fastqc_files,fastqc_files2

    """
    fastqc --nogroup -q ${in_fastq.get(0)} ${in_fastq.get(1)}
    """
}

// TODO: this may be a step we make optional, it may take a long time if the
// sample # are huge (100s or 1000s). Maybe an overall quality score heat map (per base)
// would be better?

process runDADA2QC {
    tag { "DADA2-FASTQ-QC" }
    publishDir "${params.outdir}/dada2-RawQC", mode: "copy", overwrite: true

    input:
    path("fastq/*") from dada2ReadPairsToDada2Qual.flatMap({ n -> n[1] }).collect()

    when:
    params.skip_dadaQC = false

    output:
    file '*.pdf'

    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")

    fns <- list.files("./fastq", full.names=TRUE)
    
    # this may switch to 'env' in the process at some point: 
    # https://www.nextflow.io/docs/latest/process.html?highlight=env#output-env
    # untested within R though

    pdf("qualities.pdf", onefile = TRUE)
    for (i in seq(1, length(fns), by = 4)) {
        pl <- plotQualityProfile(fns[i:(i+3)])
        print(pl)
    }

    dev.off()
    """
}

// TODO: combine MultiQC reports and split by directory (no need for two)
process runMultiQC {
    tag { "MultiQC" }
    publishDir "${params.outdir}/MultiQC-Raw", mode: 'copy', overwrite: true

    input:
    file('./raw-seq/*') from fastqc_files.collect()

    output:
    file "*_report.html" into multiqc_report
    file "*_data"

    when:
    params.skip_multiQC = false

    script:
    interactivePlots = params.interactiveMultiQC == true ? "-ip" : ""
    """
    multiqc ${interactivePlots} .
    """
}

/* ITS amplicon filtering */

// Note: should explore cutadapt options more: https://github.com/benjjneb/dada2/issues/785
// https://cutadapt.readthedocs.io/en/stable/guide.html#more-than-one

// if (params.amplicon == 'ITS') {

//     process itsFilterAndTrimStep1 {
//         tag { "ITS_Step1_${pairId}" }

//         input:
//         set pairId, file(reads) from dada2ReadPairs

//         output:
//         set val(pairId), "${pairId}.R[12].noN.fastq.gz" optional true into itsStep2
//         set val(pairId), "${pairId}.out.RDS" into itsStep3Trimming  // needed for join() later
//         file('forward_rc') into forwardP
//         file('reverse_rc') into reverseP

//         when:
//         params.precheck == false

//         script:
//         """
//         #!/usr/bin/env Rscript
//         library(dada2); packageVersion("dada2")
//         library(ShortRead); packageVersion("ShortRead")
//         library(Biostrings); packageVersion("Biostrings")

//         #Filter out reads with N's
//         out1 <- filterAndTrim(fwd = "${reads[0]}",
//                             filt = paste0("${pairId}", ".R1.noN.fastq.gz"),
//                             rev = "${reads[1]}",
//                             filt.rev = paste0("${pairId}", ".R2.noN.fastq.gz"),
//                             maxN = 0,
//                             multithread = ${task.cpus})
//         FWD.RC <- dada2:::rc("${params.fwdprimer}")
//         REV.RC <- dada2:::rc("${params.revprimer}")
        
//         # this may switch to 'env' in the process at some point: 
//         # https://www.nextflow.io/docs/latest/process.html?highlight=env#output-env
//         # untested within R though
        
//         forP <- file("forward_rc")
//         writeLines(FWD.RC, forP)
//         close(forP)

//         revP <- file("reverse_rc")
//         writeLines(REV.RC, revP)
//         close(revP)
        
//         saveRDS(out1, "${pairId}.out.RDS")
//         """
//     }
    
//     process itsFilterAndTrimStep2 {
//         tag { "ITS_Step2_${pairId}" }
//         publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

//         input:
//         set pairId, reads from itsStep2
//         file(forP) from forwardP
//         file(revP) from reverseP
        
//         output:
//         set val(pairId), "${pairId}.R[12].cutadapt.fastq.gz" optional true into itsStep3
//         file "*.cutadapt.out" into cutadaptToMultiQC

//         when:
//         params.precheck == false

//         script:
//         """
//         FWD_PRIMER=\$(<forward_rc)
//         REV_PRIMER=\$(<reverse_rc)
        
//         cutadapt -g "${params.fwdprimer}" -a \$FWD_PRIMER \\
//             -G "${params.revprimer}" -A \$REV_PRIMER \\
//             --cores ${task.cpus} \\
//             -n 2 \\
//             -o "${pairId}.R1.cutadapt.fastq.gz" \\
//             -p "${pairId}.R2.cutadapt.fastq.gz" \\
//             "${reads[0]}" "${reads[1]}" > "${pairId}.cutadapt.out"
//         """
//     }

//     process itsFilterAndTrimStep3 {
//         tag { "ITS_Step3_${pairId}" }
//         publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

//         input:
//         set pairId, file(reads), file(trimming) from itsStep3.join(itsStep3Trimming)

//         output:
//         set val(pairId), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" optional true into filteredReadsforQC, filteredReads
//         file "*.R1.filtered.fastq.gz" optional true into forReads
//         file "*.R2.filtered.fastq.gz" optional true into revReads
//         file "*.trimmed.txt" into trimTracking

//         when:
//         params.precheck == false

//         script:
//         """
//         #!/usr/bin/env Rscript
//         library(dada2); packageVersion("dada2")
//         library(ShortRead); packageVersion("ShortRead")
//         library(Biostrings); packageVersion("Biostrings")

//         out1 <- readRDS("${trimming}")
//         out2 <- filterAndTrim(fwd = paste0("${pairId}",".R1.cutadapt.fastq.gz"),
//                             filt = paste0("${pairId}", ".R1.filtered.fastq.gz"),
//                             rev = paste0("${pairId}",".R2.cutadapt.fastq.gz"),
//                             filt.rev = paste0("${pairId}", ".R2.filtered.fastq.gz"),
//                             maxEE = c(${params.maxEEFor},${params.maxEERev}),
//                             truncLen = c(${params.truncFor},${params.truncRev}),
//                             truncQ = ${params.truncQ},
//                             maxN = ${params.maxN},
//                             rm.phix = as.logical(${params.rmPhiX}),
//                             maxLen = ${params.maxLen},
//                             minLen = ${params.minLen},
//                             compress = TRUE,
//                             verbose = TRUE,
//                             multithread = ${task.cpus})
//         #Change input read counts to actual raw read counts
//         out3 <- cbind(out1, out2)
//         colnames(out3) <- c('input', 'filterN', 'cutadapt', 'filtered')
//         write.csv(out3, paste0("${pairId}", ".trimmed.txt"))
//         """
//     }
    
// }
// /* 16S amplicon filtering */
// else if (params.amplicon == '16S'){
process filterAndTrim {
    tag { "16s_${pairId}" }
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

    input:
    set pairId, file(reads) from dada2ReadPairs

    output:
    set val(pairId), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" optional true into filteredReadsforQC,filteredReads
    file "*.R1.filtered.fastq.gz" optional true into forReadsLE
    file "*.R2.filtered.fastq.gz" optional true into revReadsLE
    file "*.trimmed.txt" into trimTracking

    when:
    params.precheck == false

    script:
    phix = params.rmPhiX ? '--rmPhiX TRUE' : '--rmPhiX FALSE'
    """
    16S_FilterAndTrim.R ${phix} --id ${pairId} \\
        --fwd ${reads[0]} \\
        --rev ${reads[1]} \\
        --cpus ${task.cpus} \\
        --trimFor ${params.trimFor} \\
        --trimRev ${params.trimRev} \\
        --truncFor ${params.truncFor} \\
        --truncRev ${params.truncRev} \\
        --truncQ ${params.truncQ} \\
        --maxEEFor ${params.maxEEFor} \\
        --maxEERev ${params.maxEERev} \\
        --maxN ${params.maxN} \\
        --maxLen ${params.maxLen} \\
        --minLen ${params.minLen}
    """
}
cutadaptToMultiQC = Channel.empty()
// } else {
//     // We need to shut this down!
//     cutadaptToMultiQC = Channel.empty()
//     filteredReads = Channel.empty()
//     filteredReadsforQC = Channel.empty()
// }

process runFastQC_postfilterandtrim {
    tag { "rFQC_post_FT.${pairId}" }
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

process runMultiQC_postfilterandtrim {
    tag { "runMultiQC_postfilterandtrim" }
    publishDir "${params.outdir}/MultiQC-Post-FilterTrim", mode: 'copy', overwrite: true

    input:
    file('./raw-seq/*') from fastqc_files2.collect()
    file('./trimmed-seq/*') from fastqc_files_post.collect()
    file('./cutadapt/*') from cutadaptToMultiQC.collect()

    output:
    file "*_report.html" into multiqc_report_post
    file "*_data"

    when:
    params.precheck == false

    script:
    interactivePlots = params.interactiveMultiQC == true ? "-ip" : ""
    """
    multiqc ${interactivePlots} .
    """
}

process mergeTrimmedTable {
    tag { "mergeTrimmedTable" }
    publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

    input:
    file trimData from trimTracking.collect()

    output:
    file "all.trimmed.csv" into trimmedReadTracking

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    trimmedFiles <- list.files(path = '.', pattern = '*.trimmed.txt')
    sample.names <- sub('.trimmed.txt', '', trimmedFiles)
    trimmed <- do.call("rbind", lapply(trimmedFiles, function (x) as.data.frame(read.csv(x))))
    colnames(trimmed)[1] <- "Sequence"
    trimmed\$SampleID <- sample.names
    write.csv(trimmed, "all.trimmed.csv", row.names = FALSE)
    """
}

/*
 *
 * Step 2: Learn error rates (run on all samples)
 *
 */

// TODO: combine For and Rev process to reduce code duplication?

process LearnErrorsFor {
    tag { "LearnErrorsFor" }
    publishDir "${params.outdir}/dada2-LearnErrors", mode: "copy", overwrite: true

    input:
    file fReads from forReadsLE.collect()

    output:
    file "errorsF.RDS" into errorsFor
    file "*.pdf"

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2);
    packageVersion("dada2")
    setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})

    # File parsing
    filtFs <- list.files('.', pattern="R1.filtered.fastq.gz", full.names = TRUE)
    sample.namesF <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
    set.seed(100)

    # Learn forward error rates
    errF <- learnErrors(filtFs, multithread=${task.cpus})

    if (as.logical('${params.qualityBinning}') == TRUE ) {
        print("Running binning correction")
        errs <- t(apply(getErrors(errF), 1, function(x) { x[x < x[40]] = x[40]; return(x)} ))
        errF\$err_out <- errs
    }

    pdf("R1.err.pdf")
    plotErrors(errF, nominalQ=TRUE)
    dev.off()

    saveRDS(errF, "errorsF.RDS")
    """
}

process LearnErrorsRev {
    tag { "LearnErrorsRev" }
    publishDir "${params.outdir}/dada2-LearnErrors", mode: "copy", overwrite: true

    input:
    file rReads from revReadsLE.collect()

    output:
    file "errorsR.RDS" into errorsRev
    file "*.pdf"

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)
    packageVersion("dada2")
    setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})

    # load error profiles

    # File parsing
    filtRs <- list.files('.', pattern="R2.filtered.fastq.gz", full.names = TRUE)
    sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
    set.seed(100)

    # Learn forward error rates
    errR <- learnErrors(filtRs, multithread=${task.cpus})

    # optional NovaSeq binning error correction
    if (as.logical('${params.qualityBinning}') == TRUE) {
        print("Running binning correction")
        errs <- t(apply(getErrors(errR), 1, function(x) { x[x < x[40]] = x[40]; return(x)} ))
        errR\$err_out <- errs
    }

    pdf("R2.err.pdf")
    plotErrors(errR, nominalQ=TRUE)
    dev.off()

    saveRDS(errR, "errorsR.RDS")
    """
}

/*
 *
 * Step 3: Dereplication, Sample Inference, Merge Pairs
 *
 */

/*
 *
 * Step 4: Construct sequence table
 *
 */

// if (params.pool == "T" || params.pool == 'pseudo') {

    // In this track we run pooling (either complete or 'pseudopooling')
filteredReads.into{ filteredReadsR1;filteredReadsR2;filteredReadsMerge }

    process InferFor {
        tag { "InferFor" }
        publishDir "${params.outdir}/dada2-Derep-Pooled", mode: "copy", overwrite: true

        input:
        file filts from filteredReadsR1
            .map{ it[1] } // keep only R1 reads, toss ID and R2
            .collect( )
        file err from errorsFor

        output:
        file "all.ddF.RDS" into dadaMergeF,dadaForReadTrackingF

        when:
        params.precheck == false

        script:
        """
        dadaPooled.R --err ${err} \\
            --read 'R1' \\
            --pool ${params.pool} \\
            --cpus ${task.cpus} \\
            --dadaOpt ${params.dadaOpt}
        """
    }

    process InferRev {
        tag { "InferRev" }
        publishDir "${params.outdir}/dada2-Derep-Pooled", mode: "copy", overwrite: true

        input:
        file filts from filteredReadsR2
            .map{ it[1] } // keep only reads, toss ID
            .collect( )
        file err from errorsRev

        output:
        file "all.ddF.RDS" into dadaMergeR,dadaForReadTrackingR

        when:
        params.precheck == false

        script:
        """
        dadaPooled.R --err ${err} \\
            --read 'R2' \\
            --pool ${params.pool} \\
            --cpus ${task.cpus} \\
            --dadaOpt ${params.dadaOpt}
        """
    }

    process mergePooled {
        tag { "mergeDadaRDS" }
        publishDir "${params.outdir}/dada2-Derep-Pooled", mode: "copy", overwrite: true

        input:
        file ddF from dadaMergeF
        file ddR from dadaMergeR
        file filts from filteredReadsMerge.collect( )

        output:
        file "seqtab.RDS" into seqTable,rawSeqTableToRename
        file "all.mergers.RDS" into mergerTracking
        file "seqtab.*"

        when:
        params.precheck == false

        script:
        '''
        MergePairs.R --cpus ${task.cpus} \\
            --minOverlap ${params.minOverlap} \\
            --maxMismatch ${params.maxMismatch} \\
            --trimOverhang ${params.trimOverhang} \\
            --minMergedLen ${params.minMergedLen} \\
            --maxMergedLen ${params.maxMergedLen} \\
            --justConcatenate ${params.justConcatenate}
        '''
    }

// } else {
//     // pool = F, process per sample
//     process PerSampleInferDerepAndMerge {
//         tag { "PerSampleInferDerepAndMerge" }
//         publishDir "${params.outdir}/dada2-Derep", mode: "copy", overwrite: true

//         input:
//         set val(pairId), file(filtFor), file(filtRev) from filteredReads
//         file errFor from errorsFor
//         file errRev from errorsRev

//         output:
//         file "seqtab.RDS" into seqTable
//         file "all.mergers.RDS" into mergerTracking
//         file "all.ddF.RDS" into dadaForReadTracking
//         file "all.ddR.RDS" into dadaRevReadTracking
//         file "all.derepF.RDS" into dadaForDerep
//         file "all.derepR.RDS" into dadaRevDerep
//         file "seqtab.*"

//         when:
//         params.precheck == false

//         script:
//         """
//         #!/usr/bin/env Rscript
//         library(dada2)
//         packageVersion("dada2")
//         setDadaOpt(${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")})

//         errF <- readRDS("${errFor}")
//         errR <- readRDS("${errRev}")
//         cat("Processing:", "${pairId}", "\\n")

//         derepF <- derepFastq("${filtFor}")

//         ddF <- dada(derepF, err=errF, multithread=${task.cpus}, pool=as.logical("${params.pool}"))

//         derepR <- derepFastq("${filtRev}")
//         ddR <- dada(derepR, err=errR, multithread=${task.cpus}, pool=as.logical("${params.pool}"))

//         merger <- mergePairs(ddF, derepF, ddR, derepR,
//             returnRejects = TRUE,
//             minOverlap = ${params.minOverlap},
//             maxMismatch = ${params.maxMismatch},
//             trimOverhang = as.logical("${params.trimOverhang}"),
//             justConcatenate=as.logical("${params.justConcatenate}")
//             )

//         saveRDS(merger, paste("${pairId}", "merged", "RDS", sep="."))

//         saveRDS(ddFs, "all.ddF.RDS")
//         saveRDS(derepFs, "all.derepFs.RDS")

//         saveRDS(ddRs, "all.ddR.RDS")
//         saveRDS(derepRs, "all.derepRs.RDS")
//         """
//     }

//     process mergeDadaRDS {
//         tag { "mergeDadaRDS" }
//         publishDir "${params.outdir}/dada2-Inference", mode: "copy", overwrite: true

//         input:
//         file ddFs from dadaFor.collect()
//         file ddRs from dadaRev.collect()

//         output:
//         file "all.ddF.RDS" into dadaForReadTracking
//         file "all.ddR.RDS" into dadaRevReadTracking

//         when:
//         params.precheck == false

//         script:
//         '''
//         #!/usr/bin/env Rscript
//         library(dada2)
//         packageVersion("dada2")

//         dadaFs <- lapply(list.files(path = '.', pattern = '.ddF.RDS$'), function (x) readRDS(x))
//         names(dadaFs) <- sub('.ddF.RDS', '', list.files('.', pattern = '.ddF.RDS'))
//         dadaRs <- lapply(list.files(path = '.', pattern = '.ddR.RDS$'), function (x) readRDS(x))
//         names(dadaRs) <- sub('.ddR.RDS', '', list.files('.', pattern = '.ddR.RDS'))
//         saveRDS(dadaFs, "all.ddF.RDS")
//         saveRDS(dadaRs, "all.ddR.RDS")
//         '''
//     }

//     process SequenceTable {
//         tag { "SequenceTable" }
//         publishDir "${params.outdir}/dada2-SeqTable", mode: "copy", overwrite: true

//         input:
//         file mr from mergedReads.collect()

//         output:
//         file "seqtab.RDS" into seqTable,rawSeqTableToRename
//         file "all.mergers.RDS" into mergerTracking

//         when:
//         params.precheck == false

//         script:
//         '''
//         #!/usr/bin/env Rscript
//         library(dada2)
//         packageVersion("dada2")
        
//         mergerFiles <- list.files(path = '.', pattern = '.*.RDS$')
//         pairIds <- sub('.merged.RDS', '', mergerFiles)
//         mergers <- lapply(mergerFiles, function (x) readRDS(x))
//         names(mergers) <- pairIds
//         seqtab <- makeSequenceTable(mergers)
//         seqtab <- seqtab[,nchar(colnames(seqtab)) >= ${params.minLen}]

//         saveRDS(seqtab, "seqtab.RDS")
//         saveRDS(mergers, "all.mergers.RDS")
//         '''
//     }
// }

/*
 *
 * Step 8: Remove chimeras
 *
 */

if (!params.skipChimeraDetection) {
    process RemoveChimeras {
        tag { "RemoveChimeras" }
        publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "copy", overwrite: true

        input:
        file st from seqTable

        output:
        file "seqtab_final.RDS" into seqTableToTax,seqTableToRename

        when:
        params.precheck == false

        script:
        chimOpts = params.removeBimeraDenovoOptions != false ? ", ${params.removeBimeraDenovoOptions}" : ''
        """
        #!/usr/bin/env Rscript
        library(dada2)
        packageVersion("dada2")

        st.all <- readRDS("${st}")

        # Remove chimeras
        seqtab <- removeBimeraDenovo(
            st.all, 
            method="consensus", 
            multithread=${task.cpus}, 
            verbose=TRUE ${chimOpts} 
            )

        saveRDS(seqtab, "seqtab_final.RDS")
        """
    }
} else {
    seqTable.into {seqTableToTax;seqTableToRename}
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
                tag { "AssignTaxSpeciesRDP" }
                publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "copy", overwrite: true

                input:
                file st from seqTableToTax
                file ref from refFile
                file sp from speciesFile

                output:
                file "tax_final.RDS" into taxFinal,taxTableToTable
                file "bootstrap_final.RDS" into bootstrapFinal

                when:
                params.precheck == false

                script:
                """
                #!/usr/bin/env Rscript
                library(dada2)
                packageVersion("dada2")

                seqtab <- readRDS("${st}")

                # Assign taxonomy
                tax <- assignTaxonomy(seqtab, "${ref}",
                                        multithread=${task.cpus},
                                        tryRC = TRUE,
                                        outputBootstraps = TRUE,
                                        minBoot = ${params.minBoot},
                                        verbose = TRUE)
                boots <- tax\$boot

                tax <- addSpecies(tax\$tax, "${sp}",
                                 tryRC = TRUE,
                                 verbose = TRUE)

                rownames(tax) <- colnames(seqtab)

                # Write original data
                saveRDS(tax, "tax_final.RDS")
                saveRDS(boots, "bootstrap_final.RDS")
                """
            }

        } else {

            process AssignTaxonomyRDP {
                tag { "TaxonomyRDP" }
                publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "copy", overwrite: true

                input:
                file st from seqTableToTax
                file ref from refFile

                output:
                file "tax_final.RDS" into taxFinal,taxTableToTable
                file "bootstrap_final.RDS" into bootstrapFinal

                when:
                params.precheck == false

                script:
                taxLevels = params.taxLevels ? "c( ${params.taxLevels} )," : ''
                """
                #!/usr/bin/env Rscript
                library(dada2)
                packageVersion("dada2")

                seqtab <- readRDS("${st}")

                # Assign taxonomy
                tax <- assignTaxonomy(seqtab, "${ref}",
                                      multithread=${task.cpus},
                                      minBoot = ${params.minBoot},
                                      tryRC = TRUE,
                                      outputBootstraps = TRUE, ${taxLevels}
                                      verbose = TRUE 
                                      )

                # Write to disk
                saveRDS(tax\$tax, "tax_final.RDS")
                saveRDS(tax\$boot, "bootstrap_final.RDS")
                """
            }
        }
    } else if (params.taxassignment == 'idtaxa') {
        // Experimental!!! This assigns full taxonomy to species level, but only for
        // some databases; unknown whether this works with concat sequences.  ITS
        // doesn't seem to be currently supported
        process TaxonomyIDTAXA {
            tag { "TaxonomyIDTAXA" }
            publishDir "${params.outdir}/dada2-Chimera-Taxonomy", mode: "copy", overwrite: true

            input:
            file st from seqTableToTax
            file ref from refFile // this needs to be a database from the IDTAXA site

            output:
            file "tax_final.RDS" into taxFinal,taxTableToTable
            file "bootstrap_final.RDS" into bootstrapFinal
            file "raw_idtaxa.RDS"

            when:
            params.precheck == false

            script:
            """
            #!/usr/bin/env Rscript
            library(dada2)
            library(DECIPHER)
            packageVersion("DECIPHER")

            seqtab <- readRDS("${st}")

            # Create a DNAStringSet from the ASVs
            dna <- DNAStringSet(getSequences(seqtab))

            # load database; this should be a RData file
            load("${refFile}")

            ids <- IdTaxa(dna, trainingSet,
                strand="both",
                processors=${task.cpus},
                verbose=TRUE)
            # ranks of interest
            ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
            saveRDS(ids, 'raw_idtaxa.RDS')

            # Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
            taxid <- t(sapply(ids, function(x) {
                    m <- match(ranks, x\$rank)
                    taxa <- x\$taxon[m]
                    taxa[startsWith(taxa, "unclassified_")] <- NA
                    taxa
            }))
            colnames(taxid) <- ranks
            rownames(taxid) <- getSequences(seqtab)

            boots <- t(sapply(ids, function(x) {
                    m <- match(ranks, x\$rank)
                    bs <- x\$confidence[m]
                    bs
            }))
            colnames(boots) <- ranks
            rownames(boots) <- getSequences(seqtab)

            # Write to disk
            saveRDS(taxid, "tax_final.RDS")
            saveRDS(boots, "bootstrap_final.RDS")
            """
        }

    } else if (params.taxassignment) {
        exit 1, "Unknown taxonomic assignment method set: ${params.taxassignment}"
    } else {
        exit 1, "No taxonomic assignment method set, but reference passed"
    }
} else {
    // set tax channels to 'false', do NOT assign taxonomy
    taxFinal = Channel.empty()
    taxTableToTable = Channel.empty()
    bootstrapFinal = Channel.empty()
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
    tag { "RenameASVs" }
    publishDir "${params.outdir}/dada2-Tables", mode: "copy", overwrite: true

    input:
    file st from seqTableToRename
    file rawst from rawSeqTableToRename

    output:
    file "seqtab_final.simple.RDS" into seqTableFinalToBiom,seqTableFinalToTax,seqTableFinalTree,seqTableFinalTracking,seqTableToTable,seqtabToPhyloseq,seqtabToTaxTable
    file "asvs.${params.idType}.nochim.fna" into seqsToAln, seqsToQIIME2
    file "readmap.RDS" into readsToRenameTaxIDs // needed for remapping tax IDs
    file "asvs.${params.idType}.raw.fna"

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)
    library(ShortRead)
    library(digest)

    # read RDS w/ data
    st <- readRDS("${st}")
    st.raw <- readRDS("${rawst}")

    # get sequences
    seqs <- colnames(st)
    seqs.raw <- colnames(st.raw)

    # get IDs based on idType
    ids_study <- switch("${params.idType}", simple=paste("ASV", 1:ncol(st), sep = ""),
                                md5=sapply(colnames(st), digest, algo="md5"))
    ids_study.raw <- switch("${params.idType}", simple=paste("ASV", 1:ncol(st.raw), sep = ""),
                                md5=sapply(colnames(st.raw), digest, algo="md5"))
    
    # sub IDs
    colnames(st) <- ids_study
    colnames(st.raw) <- ids_study.raw

    # generate FASTA
    seqs.dna <- ShortRead(sread = DNAStringSet(seqs), id = BStringSet(ids_study))
    # Write out fasta file.
    writeFasta(seqs.dna, file = 'asvs.${params.idType}.nochim.fna')

    seqs.dna.raw <- ShortRead(sread = DNAStringSet(seqs.raw), id = BStringSet(ids_study.raw))
    writeFasta(seqs.dna.raw, file = 'asvs.${params.idType}.raw.fna')

    # Write modified data (note we only keep the no-chimera reads for the next stage)
    saveRDS(st, "seqtab_final.simple.RDS")
    saveRDS(data.frame(id = ids_study, seq = seqs), "readmap.RDS")
    """
}

process GenerateSeqTables {
    tag { "GenerateSeqTables" }
    publishDir "${params.outdir}/dada2-Tables", mode: "link", overwrite: true

    input:
    file st from seqTableToTable

    output:
    file "seqtab_final.simple.qiime2.txt" into featuretableToQIIME2
    file "*.txt"

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)
    library(ShortRead)

    seqtab <- readRDS("${st}")

    if (as.logical('${params.sampleRegex}' != FALSE )) {
        rownames(seqtab) <- gsub('${params.sampleRegex}', "\\\\1", rownames(seqtab), perl = TRUE)
    }

    # Generate table output
    write.table(data.frame('SampleID' = row.names(seqtab), seqtab),
        file = 'seqtab_final.txt',
        row.names = FALSE,
        col.names=c('#SampleID', colnames(seqtab)), sep = "\t")

    ######################################################################
    # Convert to simple table + FASTA, from
    # https://github.com/LangilleLab/microbiome_helper/blob/master/convert_dada2_out.R#L69
    ######################################################################

    # Generate OTU table output (rows = samples, cols = ASV)
    write.table(data.frame('SampleID' = row.names(seqtab), seqtab),
        file = 'seqtab_final.simple.txt',
        row.names = FALSE,
        col.names=c('#SampleID', colnames(seqtab)),
        sep = "\t")

    # Generate OTU table for QIIME2 import (rows = ASVs, cols = samples)
    write.table(
        data.frame('Taxa' = colnames(seqtab), t(seqtab), check.names = FALSE),
        file = 'seqtab_final.simple.qiime2.txt',
        row.names = FALSE,
        quote=FALSE,
        sep = "\t")

    # Write modified data
    saveRDS(seqtab, "seqtab_final.simple.RDS")
    """
}

process GenerateTaxTables {
    tag { "GenerateTaxTables" }
    publishDir "${params.outdir}/dada2-Tables", mode: "link", overwrite: true

    input:
    file tax from taxTableToTable
    file bt from bootstrapFinal
    file map from readsToRenameTaxIDs

    output:
    file "tax_final.simple.RDS" into taxtabToPhyloseq
    file "tax_final.simple.txt" into taxtableToQIIME2
    file "*.txt"

    when:
    params.precheck == false

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2)
    library(ShortRead)

    tax <- readRDS("${tax}")
    map <- readRDS("${map}")

    # Note that we use the old ASV ID for output here
    write.table(data.frame('ASVID' = row.names(tax), tax),
        file = 'tax_final.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(tax)), sep = "\t")

    # Tax table
    if(!identical(rownames(tax), as.character(map\$seq))){
        stop("sequences in taxa and sequence table are not ordered the same.")
    }

    tax[is.na(tax)] <- "Unclassified"
    rownames(tax) <- map\$id
    taxa_combined <- apply(tax, 1, function(x) paste(x, collapse=";"))
    taxa_out <- data.frame(names(taxa_combined), taxa_combined)
    colnames(taxa_out) <- c("#OTU ID", "taxonomy")

    write.table(data.frame('ASVID' = row.names(tax), tax),
        file = 'tax_final.simple.full.txt',
        row.names = FALSE,
        col.names=c('#OTU ID', colnames(tax)), sep = "\t")

    write.table(taxa_out,
        file = 'tax_final.simple.txt',
        row.names = FALSE,
        sep = "\t")

    if (file.exists('bootstrap_final.RDS')) {
        boots <- readRDS("${bt}")
        if(!identical(rownames(boots), as.character(map\$seq))){
            stop("sequences in bootstrap and sequence table are not ordered the same.")
        }
        rownames(boots) <- map\$id
        write.table(data.frame('ASVID' = row.names(boots), boots),
            file = 'tax_final.bootstraps.simple.full.txt',
            row.names = FALSE,
            col.names=c('#OTU ID', colnames(boots)), sep = "\t")
    }

    # Write modified data
    saveRDS(tax, "tax_final.simple.RDS")
    """
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
            tag { "AlignReadsInfernal" }
            publishDir "${params.outdir}/dada2-Infernal", mode: "copy", overwrite: true

            input:
            file seqs from seqsToAln
            file cm from cmFile

            output:
            file "aligned_seqs.stk"
            file "aln.scores"
            file "aligned_seqs.fasta" into alnFile,alnToQIIME2

            script:
            """
            # from the original IM-TORNADO pipeline
            cmalign --cpu ${task.cpus} \\
                  -g --notrunc --sub --dnaout --noprob \\
                  --sfile aln.scores \\
                  -o aligned_seqs.stk \\
                  ${cm} ${seqs}

            # script from P. Jeraldo (Mayo)
            stkToFasta.py aligned_seqs.stk aligned_seqs.fasta
            """
        }
    } else if (params.aligner == 'DECIPHER') {

        process AlignReadsDECIPHER {
            tag { "AlignReadsDECIPHER" }
            publishDir "${params.outdir}/dada2-DECIPHER", mode: "copy", overwrite: true
            errorStrategy 'ignore'

            input:
            file seqs from seqsToAln

            output:
            file "aligned_seqs.fasta" optional true into alnFile,alnToQIIME2
            
            script:
            """
            #!/usr/bin/env Rscript
            library(dada2)
            library(DECIPHER)

            seqs <- readDNAStringSet("${seqs}")
            alignment <- AlignSeqs(seqs,
                                   anchor=NA,
                                   processors = ${task.cpus})
            writeXStringSet(alignment, "aligned_seqs.fasta")
            """
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
            tag { "GenerateTreePhangorn" }
            publishDir "${params.outdir}/dada2-Phangorn", mode: "copy", overwrite: true

            input:
            file aln from alnFile

            output:
            file "phangorn.tree.RDS" into treeRDS
            file "tree.newick" into treeFile
            file "tree.GTR.newick" into treeGTRFile

            script:
            """
            #!/usr/bin/env Rscript
            library(phangorn)

            phang.align <- read.phyDat("aligned_seqs.fasta",
                                        format = "fasta",
                                        type = "DNA")

            dm <- dist.ml(phang.align)
            treeNJ <- NJ(dm) # Note, tip order != sequence order
            fit = pml(treeNJ, data=phang.align)
            write.tree(fit\$tree, file = "tree.newick")

            ## negative edges length changed to 0!
            fitGTR <- update(fit, k=4, inv=0.2)
            fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                                  rearrangement = "stochastic", control = pml.control(trace = 0))
            saveRDS(fitGTR, "phangorn.tree.RDS")
            write.tree(fitGTR\$tree, file = "tree.GTR.newick")
            """
        }
    } else if (params.runTree == 'fasttree') {

        process GenerateTreeFasttree {
            tag { "GenerateTreeFasttree" }
            publishDir "${params.outdir}/dada2-Fasttree", mode: "copy", overwrite: true

            input:
            file aln from alnFile

            output:
            file "fasttree.tree" into treeGTRFile, treeToQIIME2
            // need to deadend the other channels, they're hanging here

            script:
            """
            OMP_NUM_THREADS=${task.cpus} FastTree -nt \\
                -gtr -gamma -spr 4 -mlacc 2 -slownni \\
                -out fasttree.tree \\
                aligned_seqs.fasta
            """
        }

    } else {
        // dead-end channels generated above
    }

    process RootTree {
        tag { "RootTree" }
        publishDir "${params.outdir}/dada2-RootedTree", mode: "link"

        input:
        file tree from treeGTRFile

        output:
        file "rooted.newick" into rootedTreeFile, rootedToQIIME2
        // need to deadend the other channels, they're hanging here

        script:
        """
        #!/usr/bin/env Rscript
        library(phangorn)
        library(ape)

        tree <- read.tree(file = "${tree}")

        midtree <- midpoint(tree)

        write.tree(midtree, file = "rooted.newick")
        """
    }
} else {
    // Note these are caught downstream
    alnToQIIME2 = Channel.empty()
    treeToQIIME2 = Channel.empty()
    rootedToQIIME2 = Channel.empty()
}

// TODO: rewrite using the python BIOM tools

process BiomFile {
    tag { "BiomFile" }
    publishDir "${params.outdir}/dada2-BIOM", mode: "copy", overwrite: true

    input:
    file sTable from seqTableFinalToBiom
    file tTable from taxFinal

    output:
    file "dada2.biom" into biomFile

    when:
    params.precheck == false & params.toBIOM == true

    script:
    """
    #!/usr/bin/env Rscript
    library(biomformat)
    packageVersion("biomformat")
    seqtab <- readRDS("${sTable}")
    taxtab <- readRDS("${tTable}")
    st.biom <- make_biom(t(seqtab), observation_metadata = taxtab)
    write_biom(st.biom, "dada2.biom")
    """
}

/*
 *
 * Step 10: Track reads
 *
 */

// Broken?: needs a left-join on the initial table

// process ReadTracking {
//     tag { "ReadTracking" }
//     publishDir "${params.outdir}/dada2-ReadTracking", mode: "copy", overwrite: true

//     input:
//     file trimmedTable from trimmedReadTracking
//     file sTable from seqTableFinalTracking
//     file mergers from mergerTracking
//     file ddFs from dadaForReadTracking
//     file ddRs from dadaRevReadTracking

//     output:
//     file "all.readtracking.txt"

//     when:
//     params.precheck == false

//     script:
//     """
//     #!/usr/bin/env Rscript
//     library(dada2)
//     packageVersion("dada2")
//     library(dplyr)

//     getN <- function(x) sum(getUniques(x))

//     # the gsub here might be a bit brittle...
//     dadaFs <- as.data.frame(sapply(readRDS("${ddFs}"), getN))
//     rownames(dadaFs) <- gsub('.R1.filtered.fastq.gz', '',rownames(dadaFs))
//     colnames(dadaFs) <- c("denoisedF")
//     dadaFs\$SampleID <- rownames(dadaFs)

//     dadaRs <- as.data.frame(sapply(readRDS("${ddRs}"), getN))
//     rownames(dadaRs) <- gsub('.R2.filtered.fastq.gz', '',rownames(dadaRs))
//     colnames(dadaRs) <- c("denoisedR")
//     dadaRs\$SampleID <- rownames(dadaRs)

//     all.mergers <- readRDS("${mergers}")
//     mergers <- as.data.frame(sapply(all.mergers, function(x) sum(getUniques(x %>% filter(accept)))))
//     rownames(mergers) <- gsub('.R1.filtered.fastq.gz', '',rownames(mergers))
//     colnames(mergers) <- c("merged")
//     mergers\$SampleID <- rownames(mergers)

//     seqtab.nochim <- as.data.frame(rowSums(readRDS("${sTable}")))
//     rownames(seqtab.nochim) <- gsub('.R1.filtered.fastq.gz', '',rownames(seqtab.nochim))
//     colnames(seqtab.nochim) <- c("seqtab.nochim")
//     seqtab.nochim\$SampleID <- rownames(seqtab.nochim)

//     trimmed <- read.csv("${trimmedTable}")

//     track <- Reduce(function(...) merge(..., by = "SampleID",  all.x=TRUE),  list(trimmed, dadaFs, dadaRs, mergers, seqtab.nochim))
//     # dropped data in later steps gets converted to NA on the join
//     # these are effectively 0
//     track[is.na(track)] <- 0
    
//     write.table(track, "all.readtracking.txt", sep = "\t", row.names = FALSE)
//     """
// }

if (params.toQIIME2) {

    process toQIIME2FeatureTable {
        tag { "QIIME2-Output" }
        label 'QIIME2'
        publishDir "${params.outdir}/dada2-QIIME2", mode: "link"

        input:
        file seqtab from featuretableToQIIME2

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
            --output-path seqtab_final.simple.qza \
            --type 'FeatureTable[Frequency]'
        """
    }

    process toQIIME2TaxTable {
        tag { "QIIME2-Output" }
        label 'QIIME2'
        publishDir "${params.outdir}/dada2-QIIME2", mode: "link"

        input:
        file taxtab from taxtableToQIIME2

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
            --output-path tax_final.simple.qza \
            --type 'FeatureData[Taxonomy]'
        """
    }

    process toQIIME2Seq {
        tag { "QIIME2-Output" }
        label 'QIIME2'
        publishDir "${params.outdir}/dada2-QIIME2", mode: "link"

        input:
        file seqs from seqsToQIIME2

        output:
        file "*.qza"

        script:
        """
        qiime tools import \
            --input-path ${seqs} \
            --output-path sequences.qza \
            --type 'FeatureData[Sequence]'
        """
    }

    process toQIIME2Aln {
        tag { "QIIME2-Output" }
        label 'QIIME2'
        publishDir "${params.outdir}/dada2-QIIME2", mode: "link"

        input:
        file aln from alnToQIIME2

        when:
        alnToQIIME2 != false

        output:
        file "*.qza"

        script:
        """
        qiime tools import \
            --input-path ${aln} \
            --output-path aligned-sequences.qza \
            --type 'FeatureData[AlignedSequence]'
        """
    }

    process toQIIME2Tree {
        tag { "QIIME2-Output" }
        label 'QIIME2'
        publishDir "${params.outdir}/dada2-QIIME2", mode: "link"

        input:
        file rooted from rootedToQIIME2
        file tree from treeToQIIME2

        output:
        file "*.qza"

        script:
        """
        qiime tools import \
            --input-path ${tree} \
            --output-path unrooted-tree.qza \
            --type 'Phylogeny[Unrooted]'

        qiime tools import \
            --input-path ${rooted} \
            --output-path rooted-tree.qza \
            --type 'Phylogeny[Rooted]'
        """
    }
}

// TODO: this could eventually go into a report process to consolidate workflow info;
//       deciding between pander, knitr, or combination of the two
process sessionInfo {
        tag { "R-sessionInfo" }
        label 'sessionInfo'
        publishDir "${params.outdir}/sessionInfo", mode: "link"

        output:
        file "sessionInfo.Rmd"

        script:
        """
        #!/usr/bin/env Rscript
        pkgs <- c(
          "RCurl",
          "tidyverse",
          "pander",
          "phangorn",
          "dplyr",
          "dada2",
          "DECIPHER",
          "digest",
          "biomformat",
          "optparse"
        )
        lapply(pkgs, require, character.only = TRUE)

        sink('sessionInfo.Rmd')
        rmd <- pander(sessionInfo(), compact = FALSE)
        sink()
        """
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
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[${params.base}/16S-rDNA-dada2-pipeline] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[${params.base}/16S-rDNA-dada2-pipeline] Sent summary e-mail to $params.email (mail)"
        }
    }
}
