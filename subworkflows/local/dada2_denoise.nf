// TODO: move to a subworkflow and implement pooled vs per-sample + optional priors
include { ILLUMINA_DADA2_LEARN_ERRORS           } from '../../modules/local/illumina_learnerrors'
include { PACBIO_DADA2_LEARN_ERRORS             } from '../../modules/local/pacbio_learnerrors'
include { DADA_INFER                            } from '../../modules/local/dadainfer'
include { POOLED_SEQTABLE                       } from '../../modules/local/pooledseqtable'

// TODO: make into a general derep/denoising subworkflow; 
// call in specific subworkflows or modules
workflow DADA2_DENOISE {

    take:
    ch_trimmed_infer
 
    main:
    ch_versions = Channel.empty()
    ch_infer = Channel.empty()
    ch_merged = Channel.empty()
    
    // START: DADA2-specific
    if (params.platform == 'pacbio') {
        PACBIO_DADA2_LEARN_ERRORS (
            ch_trimmed_infer
        )
        ch_infer = PACBIO_DADA2_LEARN_ERRORS.out.error_models.join(ch_trimmed_infer)
    }  else if (params.platform == 'illumina') {
        ILLUMINA_DADA2_LEARN_ERRORS (
            ch_trimmed_infer
        )
        ch_infer = ILLUMINA_DADA2_LEARN_ERRORS.out.error_models.join(ch_trimmed_infer)
    }

    // TODO: add single-sample ('big data') run
    // this is always in pooled mode at the moment, should be adjusted

    // deal with priors here, which are optional inputs
    for_priors = params.for_priors ? file(params.for_priors, checkIfExists: true) : file("${projectDir}/assets/dummy_file")
    rev_priors = params.rev_priors ? file(params.rev_priors, checkIfExists: true) : file("${projectDir}/assets/dummy_file")

    // if (params.pool == "T" || params.pool == 'pseudo') { 
    DADA_INFER(
        ch_infer
    )

    ch_trimmed = ch_trimmed_infer
        .map { it[1] }
        .flatten()
        .collect()

    POOLED_SEQTABLE(
        DADA_INFER.out.inferred.collect(),
        ch_trimmed)

    emit:
    inferred = DADA_INFER.out.inferred
    merged_seqs = POOLED_SEQTABLE.out.merged_seqs
    filtered_seqtable = POOLED_SEQTABLE.out.filtered_seqtable
}

