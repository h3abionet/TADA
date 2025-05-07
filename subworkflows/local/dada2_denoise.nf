include { ILLUMINA_DADA2_LEARN_ERRORS           } from '../../modules/local/illumina_learnerrors'
include { PACBIO_DADA2_LEARN_ERRORS             } from '../../modules/local/pacbio_learnerrors'
include { DADA2_POOLED_DENOISE                  } from '../../subworkflows/local/dada2_pooled_denoise'
include { DADA2_PER_SAMPLE_DENOISE              } from '../../subworkflows/local/dada2_per_sample_denoise'

workflow DADA2_DENOISE {

    take:
    ch_trimmed_batch
    ch_trimmed_parallel
 
    main:
    ch_versions = Channel.empty()
    ch_errs = Channel.empty()
    ch_inferred = Channel.empty()
    ch_filtered_seqtab = Channel.empty()
    ch_merged = Channel.empty()
    
    // START: DADA2-specific
    if (params.platform == 'pacbio') {
        PACBIO_DADA2_LEARN_ERRORS (
            ch_trimmed_batch
        )
        ch_errs = PACBIO_DADA2_LEARN_ERRORS.out.error_models
    }  else if (params.platform == 'illumina') {
        ILLUMINA_DADA2_LEARN_ERRORS (
            ch_trimmed_batch
        )
        ch_errs = ILLUMINA_DADA2_LEARN_ERRORS.out.error_models
    }

    // deal with priors here, which are optional inputs
    for_priors = params.for_priors ? file(params.for_priors, checkIfExists: true) : file("${projectDir}/assets/dummy_file")
    rev_priors = params.rev_priors ? file(params.rev_priors, checkIfExists: true) : file("${projectDir}/assets/dummy_file")

    // TODO: we should disallow running "parallel-pseudo", or warn these will
    // be ignored if they are set in favor of round 1 'pseudo' priors
    // if (params.pool == "parallel" || params.pool == "parallel-pseudo") {
    if (params.pool == "parallel" || params.pool == "parallel-pseudo") {
        per_sample_errs = ch_errs.map {it -> it[1]}.collect()
        DADA2_PER_SAMPLE_DENOISE(
            per_sample_errs,
            ch_trimmed_parallel,
            for_priors,
            rev_priors)
        
        ch_merged = DADA2_PER_SAMPLE_DENOISE.out.merged_seqs
        ch_inferred = DADA2_PER_SAMPLE_DENOISE.out.inferred
        ch_filtered_seqtab = DADA2_PER_SAMPLE_DENOISE.out.filtered_seqtable
    } else {
        // this runs standard denoising using DADA2, 
        // which sequentially processes data
        //
        // TODO: can we even use priors with 'true' or 'pseudo'?
        ch_batch_errs = ch_errs.join(ch_trimmed_batch)
        DADA2_POOLED_DENOISE(
            ch_batch_errs, ch_trimmed_batch
        )

        ch_merged = DADA2_POOLED_DENOISE.out.merged_seqs
        ch_inferred = DADA2_POOLED_DENOISE.out.inferred
        ch_filtered_seqtab = DADA2_POOLED_DENOISE.out.filtered_seqtable

    }

    // ch_versions = ch_versions.mix(DADA2_POOLED_DENOISE.out.ch_versions)

    emit:
    inferred = ch_inferred
    merged_seqs = ch_merged
    filtered_seqtable = ch_filtered_seqtab
    versions = ch_versions
}

