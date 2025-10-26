include { DADA2_LEARN_ERRORS                    } from '../../modules/local/learnerrors'
include { DADA2_DEREP_SEQS                      } from '../../modules/local/dada2derepseqs'
include { DADA2_POOLED_DENOISE                  } from '../../subworkflows/local/dada2_pooled_denoise'
include { DADA2_PER_SAMPLE_DENOISE              } from '../../subworkflows/local/dada2_per_sample_denoise'

workflow DADA2_DENOISE {

    take:
    ch_trimmed
 
    main:
    ch_versions = Channel.empty()

    ch_errs = Channel.empty()
    ch_dereps = Channel.empty()
    ch_inferred = Channel.empty()
    ch_filtered_seqtab = Channel.empty()
    ch_merged = Channel.empty()
    
    // START: DADA2-specific
    DADA2_DEREP_SEQS(ch_trimmed)

    ch_dereps = DADA2_DEREP_SEQS.out.derep_rds

    // Channel setup

    // We need to group data depending on which downstream steps are needed.  There
    // are two combinations possible

    // 1. The immediate downstream QC steps can use the meta info and the read pairs.
    //    Instead of doing handstands reusing the two channels above, we emit channels 
    //    with the reads paired if needed. ch_trimmed_parallel

    // 2. LearnErrors and the pooled denoising branch requires all R1 in one run and all 
    //    R2 in another, but the two groups can be processed in parallel.  So we set up
    //    the channels with this in mind. No sample ID info is really needed.
    //    For this 'batch' run, we use two channels combining the data and 
    //    include whether they are R1 or R1 (a 'readmode') to distinguish them. ch_trimmed_batch

    ch_dereps_per_read = ch_dereps
        .map { 
            [ 'R1', it[0].single_end ? it[1] : it[1][0] ]
        }
        .concat(
            ch_dereps
                .filter { !it[0].single_end }
                .map {
                    [ 'R2', it[1][1] ]
                }
        )
        .groupTuple(sort: true)
        .dump()

    DADA2_LEARN_ERRORS(ch_dereps_per_read) 
    ch_errs = DADA2_LEARN_ERRORS.out.error_models

    // deal with priors here, which are optional inputs
    for_priors = params.for_priors ? file(params.for_priors, checkIfExists: true) : file("${projectDir}/assets/dummy_file")
    rev_priors = params.rev_priors ? file(params.rev_priors, checkIfExists: true) : file("${projectDir}/assets/dummy_file")

    // TODO: we should disallow running "parallel-pseudo", or warn these will
    // be ignored if they are set in favor of round 1 'pseudo' priors
    if (params.pool == "parallel" || params.pool == "parallel-pseudo") {
        per_sample_errs = ch_errs.map {it -> it[1]}.collect()
        if (params.pool == "parallel") {
            DADA2_PER_SAMPLE_DENOISE(
                per_sample_errs,
                ch_trimmed_parallel,
                for_priors,
                rev_priors)
            ch_merged = DADA2_PER_SAMPLE_DENOISE.out.merged_seqs
            ch_inferred = DADA2_PER_SAMPLE_DENOISE.out.inferred
            ch_filtered_seqtab = DADA2_PER_SAMPLE_DENOISE.out.filtered_seqtable
        } else {
            error "parallel-pseudo pooling not supported yet" 
        }
    } else {
        // TODO: can we even use priors with 'true' or 'pseudo'?

        // this runs standard denoising using DADA2, which sequentially processes data
        // Currently this is the only way to run pooling

        // We also stash a copy of the full RDS with *all* derep seqs
        // to make sure these are consistent between the learnErrors step
        // and this step, which optionally pools them. For really large runs
        // this will use a ton of memory

        DADA2_POOLED_DENOISE(
            ch_errs,
            DADA2_LEARN_ERRORS.out.dereps_full
        )

        ch_merged = DADA2_POOLED_DENOISE.out.merged_seqs
        ch_inferred = DADA2_POOLED_DENOISE.out.inferred
        ch_filtered_seqtab = DADA2_POOLED_DENOISE.out.filtered_seqtable
        ch_versions = ch_versions.mix(DADA2_POOLED_DENOISE.out.versions)
    }

    emit:
    inferred = ch_inferred
    merged_seqs = ch_merged
    filtered_seqtable = ch_filtered_seqtab
    versions = ch_versions
}

