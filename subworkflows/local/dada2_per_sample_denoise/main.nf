include { PER_SAMPLE_INFER             } from '../../../modules/local/persampleinferderepmerge'
include { PER_SAMPLE_MERGE             } from '../../../modules/local/persamplemergedadards'
include { PER_SAMPLE_SEQTABLE          } from '../../../modules/local/persampleseqtable'

workflow DADA2_PER_SAMPLE_DENOISE {

    take:
    ch_errs
    ch_dereps
    for_priors
    rev_priors
    stage

    main:
    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow
    PER_SAMPLE_INFER(
        ch_dereps,
        ch_errs,
        for_priors,
        rev_priors)

    ch_dds = PER_SAMPLE_INFER.out.dds
                  .map { it[1] }
                  .flatten()
                  .collect()

    PER_SAMPLE_MERGE(
        ch_dds,
        stage)

    ch_combined = PER_SAMPLE_INFER.out.combinedReads.collect()
    ch_readmode = PER_SAMPLE_INFER.out.readmode.first()

    PER_SAMPLE_SEQTABLE(
        ch_combined,
        ch_readmode,
        stage)

    emit:
    inferred = PER_SAMPLE_MERGE.out.inferred
    merged_seqs = PER_SAMPLE_SEQTABLE.out.merged_seqs
    filtered_seqtable = PER_SAMPLE_SEQTABLE.out.filtered_seqtable
    for_priors = PER_SAMPLE_MERGE.out.priors_for
    rev_priors = PER_SAMPLE_MERGE.out.priors_rev
    versions = ch_versions                     // channel: [ versions.yml ]
}

