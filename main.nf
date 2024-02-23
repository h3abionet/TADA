#!/usr/bin/env nextflow

/*
========================================================================================
               D A D A 2   P I P E L I N E
========================================================================================
 DADA2 NEXTFLOW PIPELINE FOR UCT CBIO, HPCBio

----------------------------------------------------------------------------------------
*/

repository = workflow.repository ? workflow.repository : 'TADA'

def helpMessage() {
    log.info"""
    ===================================
     ${repository}  ~  version ${params.version}
    ===================================
    Usage:

    This pipeline can be run specifying parameters in a config file or with command line flags.
    The typical example for running the pipeline with command line flags is as follows:
    
      nextflow run h3abionet/TADA --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 \
        --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex

    The typical command for running the pipeline with your own config (instead of command line flags) is as follows:
    
      nextflow run h3abionet/TADA -c dada2_user_input.config -profile uct_hex
    
    where 'dada2_user_input.config' is the configuration file (see example 'dada2_user_input.config')
    
    NB: '-profile uct_hex' still needs to be specified from the command line

    Parameters
    ----------

    Mandatory arguments:
      -profile                      Hardware config to use. Currently profile available for UCT's HPC 'uct_hex' and UIUC's 'uiuc_singularity' - create your own if necessary
                                    NB -profile should always be specified on the command line, not in the config file      
      
    Input (mandatory): Additionally, only one of the following must be specified:
      --reads                       Path to FASTQ read input data.  If the data are single-end, set '--single-end' to true.
      --input                       Path to a sample sheet (CSV); sample sheet columns must have a headers with 'id,fastq_1,fastq_2'.  
      --seqTables                   Path to input R/dada2 sequence tables. Only sequence tables with the original ASV sequences as the identifier are supported

    Output location:
      --outdir                      The output directory where the results will be saved

    Sequencing platform and primer pair:
      --platform                    One of 'illumina' or 'pacbio'. String. 
      --primer_pair                 Primer combination or kit name.  Setting this can (in some cases) set
                                    other parameters optimized for that particular primer combination or kit name.  
                                    At the moment this only handles 'pacbio' and 'shoreline'

    Read preparation parameters:
      --trimFor                     Headcrop of read1 (set 0 if no trimming is needed). integer. 
      --trimRev                     Headcrop of read2 (set 0 if no trimming is needed). integer. 
      --truncFor                    Truncate read1 here (i.e. if you want to trim 10bp off the end of a 250bp R1, truncFor should be set to 240). Enforced before trimFor/trimRev. integer. 
      --truncRev                    Truncate read2 here ((i.e. if you want to trim 10bp off the end of a 250bp R2, truncRev should be set to 240). Enforced before trimFor/trimRev. integer. 
      --maxEEFor                    After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2. integer. 
      --maxEERev                    After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2. integer. 
      --truncQ                      Truncate reads at the first instance of a quality score less than or equal to truncQ; default=2. integer. 
      --maxN                        Discard reads with more than maxN number of Ns in read; default=0. integer. 
      --maxLen                      Maximum length of trimmed sequence; maxLen is enforced before trimming and truncation; default=Inf (no maximum). integer. 
      --minLen                      Minimum length enforced after trimming and truncation; default=50. integer. 
      --rmPhiX                      Remove PhiX from read. {"T","F"}. 

      In addition due to modifications needed for variable-length sequences (ITS), the following are also supported.  Note if these are set,
      one should leave '--trimFor/--trimRev' set to 0.

      --fwdprimer                   Provided when sequence-specific trimming is required (e.g. ITS sequences using cutadapt).  Experimental
      --revprimer                   Provided when sequence-specific trimming is required (e.g. ITS sequences using cutadapt).  Experimental

    Read merging:
      --minOverlap                  integer. Minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 integer. The maximum mismatches allowed in the overlap region; default=0
      --trimOverhang                {"T","F"}. If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off.
                                    "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction                                               primer region. Default="F" (false)

    Error models:
      --qualityBinning              Binned quality correction (e.g. NovaSeq/NextSeq).  default: false
      --errorModel                  NYI. Error model to use (one of 'illumina', 'illumina-binned', 'pacbio-ccs', 'custom'). This will replace
                                    '--qualityBinning'

    Denoising using dada:
      --dadaOpt.XXX                 Set as e.g. --dadaOpt.HOMOPOLYMER_GAP_PENALTY=-1 Global defaults for the dada function, see ?setDadaOpt in R for available options and their defaults
      --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are
                                    pseudo pooling: "pseudo", true: "T", false: "F"

    Merging arguments (optional):
      --minOverlap                  The minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 The maximum mismatches allowed in the overlap region; default=0.
      --trimOverhang                If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off. "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen
                                    when reads are longer than the amplicon and read into the other-direction primer region. Default="F" (false)
      --minMergedLen                Minimum length of fragment *after* merging; default = 0 (no minimum)
      --maxMergedLen                Maximum length of fragment *after* merging; default = 0 (no maximum)

    ASV identifiers:
      --idType                      The ASV IDs are renamed to simplify downstream analysis, in particular with downstream tools.  The
                                    default is "md5" which will run MD5 on the sequence and generate a QIIME2-like unique hash.  Alternatively, 
                                    this can be set to "ASV", which simply renames the sequences in sequencial order.  

    Taxonomic arguments.  If unset, taxonomic assignment is skipped
      --taxassignment               Taxonomic assignment method.  default = 'rdp'
      --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz). default = false
      --species                     Specify path to fasta file. See dada2 addSpecies() for more detail. default = false
      --minBoot                     Minimum bootstrap value.  default = 50
      --taxLevels                   Listing of taxonomic levels for 'assignTaxonomy'. Experimental.

    Chimera detection:
      --skipChimeraDetection        Skip chimera detection/removal; default = false
      --removeBimeraDenovoOpts      Additional removeBimeraDenovo options; default = ''

    ASV multiple sequence alignment:
      --skipAlignment               Skip alignment step; note this also skips ML phylogenetic analysis. default = false
      --aligner                     Aligner to use, options are 'DECIPHER' or 'infernal'. default = 'DECIPHER'
      --infernalCM                  Covariance model (Rfam-compliant) to use.  default = false.

    Phylogenetic analysis:
      --runTree                     Tool for ML phylogenetic analysis.  Options are 'phangorn' and 'fasttree'. default = 'phangorn'

    Additional output:
      --toBIOM                      Generate a BIOM v1 compliant output. default = true
      --toQIIME2                    Generate QZA artifacts for all data for use in QIIME2. default = false

    Sample names:
      --sampleRegex                 Modify sample names based on a regular expression. default = false.  Note this option
                                    is deprecated in favor of using a sample sheet.

    Additional options:
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run
                                    sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    Help:
      --help                        Will print out summary above when executing nextflow run uct-cbio/16S-rDNA-dada2-pipeline

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
    custom_runName = workflow.runName
}

log.info "==================================="
log.info " ${params.base}/TADA  ~  version ${params.version}"
log.info "==================================="
def summary = [:]
summary['project']        = params.project
summary['clusterOptions'] = params.clusterOptions
summary['name']           = params.name
summary['base']           = params.base
summary['version']        = params.version //pipeline version
summary['runName']        = custom_runName ?: workflow.runName
summary['reads']          = params.reads
summary['input']          = params.input
summary['seqTables']      = params.seqTables
summary['single_end']     = params.single_end
summary['outdir']         = params.outdir
summary['skip_dadaQC']    = params.skip_dadaQC  // set to run this step by default, this can fail with large sample #'s
summary['skip_multiQC']   = params.skip_multiQC  // set to run this step by default, this can fail with large sample #'s
summary['precheck']       = params.precheck
summary['amplicon']       = params.amplicon // Deprecated, to be replaced by more general settings
summary['fwdprimer']      = params.fwdprimer
summary['revprimer']      = params.revprimer
summary['skip_trimming']  = params.skip_trimming // NYI; this bypasses all trimming and QC, assumes primers are removed and sequence(s) ready for DADA2
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
summary['dadaParams']     = params.dadaParams // Check this!!! if set, these are additional arguments passed to the dada() function in the PacBio workflow 
summary['dadaOpt']        = params.dadaOpt // Check this!!!
summary['pool']           = params.pool
summary['qualityBinning'] = params.qualityBinning  // false, set to true if using binned qualities (NovaSeq)
summary['minOverlap']     = params.minOverlap
summary['maxMismatch']    = params.maxMismatch
summary['trimOverhang']   = params.trimOverhang
summary['maxMergedLen']   = params.maxMergedLen // Only run if set > 1
summary['minMergedLen']   = params.minMergedLen // Only run if set > 1
summary['justConcatenate'] = params.justConcatenate  // TODO: test using false instead of string
summary['rescueUnmerged'] = params.rescueUnmerged
summary['skipChimeraDetection'] = params.skipChimeraDetection
summary['removeBimeraDenovoOptions'] = params.removeBimeraDenovoOptions
summary['reference']      = params.reference
summary['species']        = params.species
summary['taxassignment']  = params.taxassignment // default: RDP classifier implementation in dada2
summary['minBoot']        = params.minBoot // default for dada2
summary['taxLevels']      = params.taxLevels
summary['skipAlignment']  = params.skipAlignment
summary['aligner']        = params.aligner // default
summary['infernalCM']     = params.infernalCM // NYI
summary['runTree']        = params.runTree // default, current alternative is 'fasttree'
summary['idType']         = params.idType
summary['sampleRegex']    = params.sampleRegex // Deprecated; use sample sheets and --input
summary['interactiveMultiQC'] = params.interactiveMultiQC
summary['toBIOM']         = params.toBIOM  // generate BIOM v1 output
summary['toQIIME2']       = params.toQIIME2  // generate QZA artifacts for QIIME2
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Working dir']    = workflow.workDir
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

// TODO: we need to validate/sanity-check more of the parameters
// Validate inputs

/*

Reworking this. 

We need to check the sequencing platform to determine how the samples are 
processed.  Supported platforms are 'illumina' and 'pacbio'

*/

def platform = ''
platform = params.platform.toLowerCase()

if (!(["illumina","pacbio","pacbio-kinnex"].contains(platform))) {
    exit 1, "Only supported platforms (--platform argument) are currently 'pacbio', 'pacbio-kinnex', or 'illumina'"
}

// ${deity} there has to be a better way to check this!
if ( (params.seqTables && params.reads) || 
     (params.input  && params.reads) || 
     (params.seqTables && params.input)) {
    exit 1, "Only one of --reads, --input, or --seqTables is allowed!"
}

// TODO: logic is a bit convoluted here....
if ( !params.skip_trimming && (!(params.fwdprimer) || !(params.revprimer)) ) {
    log.info "Both --fwdprimer and --revprimer should be set unless skipping all trimming steps.\n" 
    log.info "[These options will become requirements in future releases]"
} else { 
    // this is a check if the primers are supplied but trimming is also set, which is almost certainly a mistake
    if ( params.trimFor != 0 || params.trimFor != 0 ) {
        log.info "trimFor and/or trimRev are set along with one or more primers (fwdprimer/revprimer).\n" 
        log.info "This will trim additional bases *in addition to* the supplied primers!"
    }
}

if (params.aligner == 'infernal' && params.infernalCM == false){
    exit 1, "Must set covariance model using --infernalCM when using Infernal"
}

if (!(['simple','md5'].contains(params.idType))) {
    exit 1, "--idType can currently only be set to 'simple' or 'md5', got ${params.idType}"
}

// Read-specific checks
if (params.reads != false) {
    Channel
    .fromFilePairs( params.reads, size: params.single_end ? 1 : 2  )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { row -> create_fastq_channel_reads(row, params.single_end) }
    .into { dada2ReadPairsToQual; dada2ReadPairsToDada2Qual; dada2ReadPairs; dada2ReadPairsToDada2QC }
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
    samples = file(params.input)

    process Check_SampleSheet {
        tag { "Check_SampleSheet" }
        publishDir "${params.outdir}/SampleSheet", mode: "copy", overwrite: true

        input:
        file(samplesheet) from samples

        when:params.input

        output:
        file('final_samplesheet.csv') into samplesheet

        """
        check_samplesheet.py ${samplesheet} final_samplesheet.csv
        """
    }

    samplesheet
        .splitCsv(header:true, sep:',')
        .map{ row -> create_fastq_channel(row) } 
        .into { dada2ReadPairsToQual; dada2ReadPairsToDada2Qual; dada2ReadPairs; dada2ReadPairsToDada2QC }
} else {
    exit 1, "Must set either --input, --reads, or --seqTables as input"
}

/*
 *
 * Step 1: Filter and trim (run per sample?)
 *
 */

if (params.reads != false || params.input != false ) { // TODO maybe we should check the channel here

    process RunFastQC {
        tag { "FastQC-${meta.id}" }
        publishDir "${params.outdir}/FastQC-Raw", mode: "copy", overwrite: true

        input:
        tuple val(meta), file(reads) from dada2ReadPairsToQual

        when:
        params.precheck | !(params.skip_FASTQC)

        output:
        file '*_fastqc.{zip,html}' into fastqc_files

        script:
        nano = params.platform == "pacbio" ? '--nano' : ''
        nogroup = params.platform == "illumina" ? '--nogroup' : ''
        """
        fastqc -t ${task.cpus} ${nano} ${nogroup} -q ${reads}
        """
    }

    // dada2ReadPairsToDada2Qual
    //   .flatMap { n -> n[1] }
    //   .map { 
    //     ['R1', n[0]], ['R2', n[1]]
    //   }
    //   .groupTuple()
    //   .set { dada2ReadPairsToDada2QC }

    process RunDADA2QC {
        tag { "RunDADA2QC" }
        publishDir "${params.outdir}/dada2-DadaQC", mode: "copy", overwrite: true

        input:
        tuple val(meta), file(reads) from dada2ReadPairsToDada2QC

        when:
        params.precheck | !(params.skip_dadaQC)

        output:
        file '*.pdf'
        file '*.RDS'

        script:
        template "DadaQC.R"
    }

    /* ITS and PacBio amplicon filtering */

    // Note: should explore cutadapt options more: https://github.com/benjjneb/dada2/issues/785
    // https://cutadapt.readthedocs.io/en/stable/guide.html#more-than-one

    if (platform == 'pacbio' || platform == 'pacbio-kinnex') {

        process PacBioTrim {
            tag { "PacBioTrim_${meta.id}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            // TODO: Note the channel name here should probably be changed
            tuple val(meta), file(reads) from dada2ReadPairs

            output:
            // tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1
            tuple val(meta), file("${meta.id}.noprimer.fastq.gz") optional true into filteredReadsToDADA2
            file("*.cutadapt.out") into cutadaptToMultiQC
            file("${meta.id}.untrimmed.fastq.gz")

            when:
            !(params.precheck)

            script:
            strictness = params.pacbio_strict_match ? '-g' : '-a'
            """
            # Logic: we should trim out the HiFi reads and require *both* primers be present (-g).
            # This should also reorient the sequence to match the primers (--rc).
            # Keep anything longer than 50bp, and allow users to filter their data by length later
            revprimer_rc=\$( echo -n ${params.revprimer} | tr "[ATGCUNYRSWKMBDHV]" "[TACGANRYSWMKVHDB]" | rev )

            cutadapt --rc \\
                ${strictness} "${params.fwdprimer}...\${revprimer_rc}" \\
                -m 50 \\
                -j ${task.cpus} \\
                --untrimmed-output "${meta.id}.untrimmed.fastq.gz" \\
                -o "${meta.id}.noprimer.fastq.gz" \\
                ${reads} > "${meta.id}.noprimer.cutadapt.out"
            """
        }

        process PacBioFilter {
            tag { "PacBioFilter_${meta.id}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            // TODO: Note the channel name here should probably be changed
            tuple val(meta), file(reads) from filteredReadsToDADA2

            output:
            // tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1
            tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1,readsToFastQC,readsToPerSample
            file "*.trimmed.txt" into trimTracking

            when:
            !(params.precheck)

            script:
            template "PacBioFilterAndTrim.R"
            
        }

        filteredReadsR2 = Channel.empty()

    } else if (platform == 'illumina' && ( params.amplicon == 'ITS' || params.amplicon == 'variable'))  {

        // TODO: note this path is only needed when using variable length sequences
        process ITSFilterAndTrimStep1 {
            tag { "ITS_Step1_${meta.id}" }

            input:
            tuple val(meta), file(reads) from dada2ReadPairs

            output:
            tuple val(meta), file("${meta.id}.R[12].noN.fastq.gz") optional true into itsStep2
            tuple val(meta), file("${meta.id}.out.RDS") into itsStep3Trimming  // needed for join() later
            file('forward_rc') into forwardP
            // TODO make this optional if data are SE
            file('reverse_rc') into reverseP

            when:
            !(params.precheck)

            script:
            template "ITSFilterAndTrimStep1.R"
        }
        
        process ITSFilterAndTrimStep2 {
            tag { "ITS_Step2_${meta.id}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            tuple(meta), file(reads) from itsStep2
            file(forP) from forwardP
            file(revP) from reverseP
            
            output:
            tuple val(meta), file("${meta.id}.R[12].cutadapt.fastq.gz") optional true into itsStep3
            file("*.cutadapt.out") into cutadaptToMultiQC

            when:
            !(params.precheck)

            script:
            outr2 = meta.single_end ? '' : "-p ${meta.id}.R2.cutadapt.fastq.gz"
            p2 = meta.single_end ? '' : "-G ${params.revprimer} -A \$REV_PRIMER"
            """
            FWD_PRIMER=\$(<forward_rc)
            REV_PRIMER=\$(<reverse_rc)
            
            cutadapt -g ${params.fwdprimer} -a \$FWD_PRIMER ${p2} \\
                --cores ${task.cpus} \\
                -n 2 \\
                -o ${meta.id}.R1.cutadapt.fastq.gz ${outr2} \\
                ${reads} > ${meta.id}.cutadapt.out
            """
        }

        process ITSFilterAndTrimStep3 {
            tag { "ITS_Step3_${meta.id}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            tuple val(meta), file(reads), file(trimming) from itsStep3.join(itsStep3Trimming)

            output:
            tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1
            tuple val(meta), file("${meta.id}.R2.filtered.fastq.gz") optional true into filteredReadsR2
            tuple val(meta), file("${meta.id}.R[12].filtered.fastq.gz") optional true into readsToFastQC,readsToPerSample
            file "*.trimmed.txt" into trimTracking

            when:
            !(params.precheck)

            script:
            template "ITSFilterAndTrimStep3.R"
        }
    }
    /* 16S amplicon filtering */
    else if (platform == 'illumina' && params.amplicon == '16S'){
        process FilterAndTrim {
            tag { "FilterAndTrim_${meta.id}" }
            publishDir "${params.outdir}/dada2-FilterAndTrim", mode: "copy", overwrite: true

            input:
            tuple val(meta), file(reads) from dada2ReadPairs

            output:
            // BIG TODO: this should be simplified into tuples with all information
            tuple val(meta), file("${meta.id}.R1.filtered.fastq.gz") optional true into filteredReadsR1
            tuple val(meta), file("${meta.id}.R2.filtered.fastq.gz") optional true into filteredReadsR2
            tuple val(meta), file("${meta.id}.R[12].filtered.fastq.gz") optional true into readsToFastQC,readsToPerSample
            // tuple val("R2"), file("${pairId}.R2.filtered.fastq.gz") optional true into revReadsLE
            file "*.trimmed.txt" into trimTracking

            when:
            !(params.precheck)

            script:
            phix = params.rmPhiX ? '--rmPhiX TRUE' : '--rmPhiX FALSE'

            template "FilterAndTrim.R"
        }
        cutadaptToMultiQC = Channel.empty()
    } else {
        // We need to shut this down!
        Channel.empty().into {cutadaptToMultiQC;filteredReads;filteredReadsforQC}
    }

    // Channel setup

    // We need to group data depending on which downstream steps are needed.  There
    // are two combinations possible

    // 1. The immediate downstream QC steps can use the meta info and the read pairs.
    //    Instead of doing handstands reusing the two channels above, we emit channels 
    //    with the reads paired if needed.

    // 2. LearnErrors and the pooled denoising branch requires all R1 and all R2, but 
    //    the two groups can be processed in parallel.  So we set up the channels with 
    //    this in mind. No sample ID info is really needed.
    filteredReadsR1
            // .concat(filteredReadsR2.ifEmpty([]))
            .map { [ 'R1', it[1]] }
            .concat(filteredReadsR2.map {['R2', it[1]] } )
            .groupTuple(sort: true)
            .into { ReadsLE; ReadsInfer; ReadsMerge }

    process RunFastQC_postfilterandtrim {
        tag { "FastQC_post_FT.${meta.id}" }
        publishDir "${params.outdir}/FastQC-Post-FilterTrim", mode: "copy", overwrite: true

        input:
        set val(meta), file(reads) from readsToFastQC

        output:
        file '*_fastqc.{zip,html}' into fastqc_files_post

        when:
        !(params.skip_FASTQC)

        script:
        nano = params.platform == "pacbio" ? '--nano' : ''
        nogroup = params.platform == "illumina" ? '--nogroup' : ''
        """
        fastqc -t ${task.cpus} ${nano} ${nogroup} -q ${reads}
        """
    }

    // TODO: change process name
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
        !(params.skip_FASTQC)

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

    process LearnErrors {
        tag { "LearnErrors:${readmode}" }
        publishDir "${params.outdir}/dada2-LearnErrors", mode: "copy", overwrite: true

        input:
        tuple val(readmode), file(reads) from ReadsLE

        output:
        // we could have a single channel here based on the pooling strategy, but
        // for now we keep these both explicitly defined
        tuple val(readmode), file("errors.${readmode}.RDS") into errorModelsPooled
        file("errors.R[12].RDS") into errorModelsPerSample
        file("${readmode}.err.pdf")

        script:
        derepreads = 100000
        dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
        if (platform == 'pacbio') {
          template "PacBioLearnErrors.R"
        } else if (platform == 'pacbio-kinnex') {
          template "PacBioKinnexLearnErrors.R"
        } else {
          template "LearnErrors.R"
        }
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
            // DADA2 step runs on all R1 and/or on all R2
            tuple val(readmode), file(err), file(reads) from errorModelsPooled
                .join(ReadsInfer)

            output:
            // Note that the mode ('merged', 'R1', 'R2') can now potentially allow SE read analysis
            file("all.dd.${readmode}.RDS") into dadaMerge,dadaToReadTracking

            when:
            params.precheck == false

            script:
            dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
            template "DadaPooled.R"
        }

        // This one is a little tricky. We can't know a priori how many instances of reads (R1 and R2) 
        // are present outside the process, but we can determine this internally within the process 
        // when we collect all of them.
        // So, here we check the size of the collected channel containing the denoised models; if 
        // there are two then this is a paired-end run, otherwise it's single-end. Logic is in the R script
        
        process PooledSeqTable {
            tag { "PooledSeqTable:${readmode}" }
            publishDir "${params.outdir}/dada2-OriginalSeqTable", mode: "copy", overwrite: true

            input:
            // we don't care about the mode here, so only get the dds (dada-inferred) RDS files
            file(dds) from dadaMerge.collect()
            // we don't care about the mode here, we only grab the reads
            file(filts) from ReadsMerge
                .map { it[1] }
                .flatten()
                .collect()

            output:
            tuple val(readmode), file("seqtab.${readmode}.RDS") into seqTable,rawSeqTableToRename
            file "all.merged.RDS" optional true into mergerTracking,mergerQC
            file "seqtab.original.*.RDS" into seqtabQC// we keep this for comparison and possible QC

            when:
            params.precheck == false

            script:
            // We could switch this to 'paired' vs 'single-end' as well
            readmode = dds.size() == 2 ? 'merged' : 'R1'
            template "SeqTables.R"
        }

    } else {

        process PerSampleInferDerepAndMerge {
            tag { "PerSampleInferDerepAndMerge:${meta.id}" }
            publishDir "${params.outdir}/dada2-Derep-Single/Per-Sample", mode: "copy", overwrite: true

            input:
            tuple val(meta), file(reads) from readsToPerSample
            file(errs) from errorModelsPerSample.collect()

            output:
            file("${meta.id}.{R1,merged}.RDS") into combinedReads
            tuple val(meta), file("${meta.id}.dd.R{1,2}.RDS") into perSampleDadaToMerge
            val(readmode) into modeSeqTable

            when:
            params.precheck == false

            script:
            dadaOpt = !params.dadaOpt.isEmpty() ? "'${params.dadaOpt.collect{k,v->"$k=$v"}.join(", ")}'" : 'NA'
            readmode = errs.size() == 2 ? 'merged' : 'R1'
            template "PerSampleDadaInfer.R"
        }

        process MergeDadaRDS {
            tag { "mergeDadaRDS" }
            publishDir "${params.outdir}/dada2-Derep-Single", mode: "copy", overwrite: true

            input:
            file(dds) from perSampleDadaToMerge
                          .map { it[1] }
                          .flatten()
                          .collect()

            output:
            file("all.dd.R{1,2}.RDS") into dadaToReadTracking

            when:
            params.precheck == false

            script:
            template "MergePerSampleDada.R"
        }

        process SequenceTable {
            tag { "SequenceTable:${readmode}" }
            publishDir "${params.outdir}/dada2-Derep-Single", mode: "copy", overwrite: true

            input:
            file(mr) from combinedReads.collect()
            val(readmode) from modeSeqTable.first()

            output:
            tuple val(readmode), file("seqtab.${readmode}.RDS") into seqTable,rawSeqTableToRename
            file "all.merged.RDS" optional true into mergerTracking,mergerQC
            file "seqtab.original.${readmode}.RDS" into seqtabQC // we keep this for comparison and possible QC
            
            when:
            params.precheck == false

            script:
            template "PerSampleSeqTable.R"
        }
    }
} else if (params.seqTables) { // TODO maybe we should check the channel here
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
    Channel.empty().into { SEChimera;RawSEChimeraToRename;trimmedReadTracking;dadaToReadTracking;mergerTracking;mergerQC }
}

/*
 *
 * Step 8: Remove chimeras
 *
 */

if (!params.skipChimeraDetection) {
    process RemoveChimeras {
        tag { "RemoveChimeras:${seqtype}" }
        publishDir "${params.outdir}/dada2-Chimera", mode: 'copy'

        input:
        tuple val(seqtype), file(st) from seqTable

        output:
        tuple val(seqtype), file("seqtab_final.${seqtype}.RDS") into seqTableToRename

        when:
        params.precheck == false

        script:
        chimOpts = params.removeBimeraDenovoOptions != false ? ", ${params.removeBimeraDenovoOptions}" : ''
        template "RemoveChimeras.R"
    }
} else {
    if (params.skipChimeraDetection == true) {
        seqTable.into {seqTableToRename}
    } else if (params.skipChimeraDetection == 'all') { // stop completely!
        // this should effectively halt the pipeline from going further
        Channel.empty()
            .into {seqTableToTax;seqTableToRename}
    }
}

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
    publishDir "${params.outdir}/dada2-Tables", mode: "copy", overwrite: true

    input:
    // This is a tricky join; we have two incoming channels that we want to
    // join by sequence type (merged, R1, R2). seqTableToRename can have all
    // three, but rawSeqTableToRename (originating from merging) has only
    // one, while SE-specific seqtables (with RawSEChimeraToRename) are 
    // generated on demand later, primarily to act as a switch.  So we join() the 
    // channels by ID, but we concatenate the original seqtable channels together 

    tuple val(seqtype), file(st), file(rawst) from seqTableToRename
        .join(rawSeqTableToRename)

    output:
    tuple val(seqtype), file("seqtab_final.${params.idType}.${seqtype}.RDS") into seqTableFinalToBiom,seqTableFinalToTax,seqTableFinalTree,seqTableFinalTracking,seqTableToTable,seqtabToPhyloseq,seqtabToTaxTable
    tuple val(seqtype), file("asvs.${params.idType}.nochim.${seqtype}.fna") into seqsToAln,seqsToQIIME2,seqsToShoreline
    tuple val(seqtype), file("readmap.${seqtype}.RDS") into seqsToTax,readsToRenameTaxIDs // needed for remapping tax IDs
    file "asvs.${params.idType}.raw.${seqtype}.fna"

    script:
    template "RenameASVs.R"
}

/*
 *
 * Step 9: Taxonomic assignment
 *
 */

if (platform == 'pacbio' && params.amplicon == 'strainid' && !params.pacbio_strict_match) {
    // Trim down StrainID ASVs to just the 16S region for taxonomic assignment.
    // This is a very specific processing step, standard 16S primers V1-V9 primers
    // are used: 
    // Shoreline 16S = For:AGRRTTYGATYHTDGYTYAG, Rev:YCNTTCCYTYDYRGTACT(rc)
    // Std = For:AGRGTTYGATYMTGGCTCAG,  Rev:RGYTACCTTGTTACGACTT(not rc!)
    process TrimStrainID {
        tag { "TrimStrainID" }
        publishDir "${params.outdir}/dada2-TrimStrainID", mode: "copy", overwrite: true

        input:
        tuple val(seqtype), file(seqs) from seqsToShoreline

        output:
        tuple val(seqtype), file("asvs.${params.idType}.nochim.${seqtype}.16S-only.fna")
        file("asvs.${params.idType}.nochim.${seqtype}.untrimmed.fna")
        file("asvs.shoreline.cutadapt.log")

        when:
        params.precheck == false

        script:
        """
        FOR=AGRGTTYGATYMTGGCTCAG
        REV=RGYTACCTTGTTACGACTT

        revprimer_rc=\$( echo -n \$REV | tr "[ATGCUNYRSWKMBDHV]" "[TACGANRYSWMKVHDB]" | rev )
        cutadapt -a \$revprimer_rc \\
                 -m 50 -j ${task.cpus} \\
                 --untrimmed-output "asvs.${params.idType}.nochim.${seqtype}.untrimmed.fna" \\
                 -o "asvs.${params.idType}.nochim.${seqtype}.16S-only.fna" \\
                 "${seqs}" > "asvs.shoreline.cutadapt.log"
        """
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
                publishDir "${params.outdir}/dada2-Taxonomy", mode: "copy"

                input:
                tuple val(seqtype), file(st) from seqsToTax
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
                publishDir path: "${params.outdir}/dada2-Taxonomy", mode: "copy", overwrite: true

                input:
                tuple val(seqtype), file(st) from seqsToTax
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
            publishDir "${params.outdir}/dada2-Taxonomy", mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(st) from seqsToTax
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

process GenerateSeqTables {
    tag { "GenerateSeqTables:${seqtype}" }
    publishDir "${params.outdir}/dada2-Tables", mode: "copy", overwrite: true

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
    publishDir "${params.outdir}/dada2-Tables", mode: "copy", overwrite: true

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
            publishDir "${params.outdir}/dada2-Infernal", mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(seqs) from seqsToAln
            file cm from cmFile

            output:
            // file "aligned_seqs.${seqtype}.stk" into alnToFASTA
            file "aln.${seqtype}.scores"

            script:
            """
            # from the original IM-TORNADO pipeline
            cmalign --cpu ${task.cpus} \\
                  -g --notrunc --sub --dnaout --noprob \\
                  --sfile aln.${seqtype}.scores \\
                  -o aligned_seqs.${seqtype}.stk \\
                  ${cm} ${seqs}
            """
        }

        process StockholmToFASTA {
            tag { "StockholmToFASTA:${seqtype}" }
            publishDir "${params.outdir}/dada2-Infernal", mode: "copy", overwrite: true

            input:
            tuple val(seqtype), file(seqs) from seqsToAln
            file cm from cmFile

            output:
            tuple val(seqtype), file("aligned_seqs.${seqtype}.fasta") optional true into alnFile,alnToQIIME2

            script:
            """
            # script from P. Jeraldo (Mayo)
            stkToFasta.py aligned_seqs.${seqtype}.stk aligned_seqs.${seqtype}.fasta
            """
        }
    } else if (params.aligner == 'DECIPHER') {

        process AlignReadsDECIPHER {
            tag { "AlignReadsDECIPHER:${seqtype}" }
            publishDir "${params.outdir}/dada2-DECIPHER", mode: "copy", overwrite: true

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
      * Step 10b: Construct phylogenetic tree
      */
     

    if (params.runTree == 'phangorn') {

        process GenerateTreePhangorn {
            tag { "GenerateTreePhangorn:${seqtype}" }
            publishDir "${params.outdir}/dada2-Phangorn", mode: "copy", overwrite: true

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
            publishDir "${params.outdir}/dada2-Fasttree", mode: "copy", overwrite: true

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
        publishDir "${params.outdir}/dada2-RootedTree", mode: "copy"

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
    publishDir "${params.outdir}/dada2-BIOM", mode: "copy", overwrite: true

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
    file(trimmedTable) from trimmedReadTracking
    tuple(seqtype), file(sTable) from seqTableFinalTracking
    file(dds) from dadaToReadTracking.collect()
    file(mergers) from mergerTracking.ifEmpty([])

    output:
    file "all.readtracking.${seqtype}.txt"

    script:
    template "ReadTracking.R"
}

process MergingQC {
    tag { "MergingQC" }
    publishDir "${params.outdir}/dada2-DadaQC", mode: "copy", overwrite: true

    input:
    file(mergers) from mergerQC

    output:
    file "read-overlap-heatmap*"

    script:
    template "PlotMergers.R"
}

process SeqLengthQC {
    tag { "SeqLengthQC" }
    publishDir "${params.outdir}/dada2-DadaQC", mode: "copy", overwrite: true

    input:
    file(seqtab) from seqtabQC

    output:
    file "asv-length-distribution*"

    script:
    template "PlotSeqLengths.R"
}

process ToQIIME2FeatureTable {
    tag { "QIIME2-SeqTable:${seqtype}" }
    label 'QIIME2'
    publishDir "${params.outdir}/dada2-QIIME2", mode: "copy"

    when: params.toQIIME2

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
    publishDir "${params.outdir}/dada2-QIIME2", mode: "copy"


    input:
    tuple val(seqtype), file(taxtab) from taxtableToQIIME2

    output:
    file "*.qza"

    when: params.toQIIME2 && taxtableToQIIME2 != false

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
    publishDir "${params.outdir}/dada2-QIIME2", mode: "copy"
    
    when: params.toQIIME2

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
    publishDir "${params.outdir}/dada2-QIIME2", mode: "copy"

    when: params.toQIIME2

    input:
    tuple val(seqtype), file(aln) from alnToQIIME2

    when: params.toQIIME2 && alnToQIIME2 != false

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
    publishDir "${params.outdir}/dada2-QIIME2", mode: "copy"

    when: params.toQIIME2

    input:
    tuple val(seqtype), file(rooted), file(tree) from rootedToQIIME2.join(treeToQIIME2)

    output:
    file "*.qza"

    script:
    """
    qiime tools import \
        --input-path ${tree} \
        --output-path unrooted-tree.${seqtype}.qza \
        --type 'Phylogeny[Unrooted]'

    qiime tools import \
        --input-path ${rooted} \
        --output-path rooted-tree.${seqtype}.qza \
        --type 'Phylogeny[Rooted]'
    """
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

    def subject = "[${params.base}/TADA] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[${params.base}/TADA] FAILED: $workflow.runName"
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
          log.info "[${params.base}/TADA] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[${params.base}/TADA] Sent summary e-mail to $params.email (mail)"
        }
    }
}

// code modified from the nf-core RNA-Seq workflow
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}

def create_fastq_channel_reads(ArrayList row, Boolean single_end=FALSE) {
    // create meta map
    def meta = [:]
    meta.id         = row[0]
    meta.single_end = single_end

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row[1][0]).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row[1][0]}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row[1][0]) ] ]
    } else {
        if (!file(row[1][1]).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row[1][1]}"
        }
        fastq_meta = [ meta, [ file(row[1][0]), file(row[1][1]) ] ]
    }
    return fastq_meta
}
