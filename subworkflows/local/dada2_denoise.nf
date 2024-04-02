// TODO: move to a subworkflow and implement pooled vs per-sample + optional priors
include { LEARN_ERRORS           } from '../../modules/local/learnerrors'
include { DADA_INFER             } from '../../modules/local/dadainfer'
include { POOLED_SEQTABLE        } from '../../modules/local/pooledseqtable'
include { REMOVE_CHIMERAS        } from '../../modules/local/removechimeras'
include { RENAME_ASVS            } from '../../modules/local/renameasvs'

workflow DADA2_DENOISE {

    take:
    // TODO nf-core: edit input (take) channels
    ch_trimmed_infer // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    LEARN_ERRORS (
        ch_trimmed_infer
    )

    ch_infer = LEARN_ERRORS.out.error_models.join(ch_trimmed_infer)

    // TODO: add single-sample ('big data') run
    // this is always in pooled mode at the moment, should be adjusted
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
        ch_trimmed
        )

    REMOVE_CHIMERAS(
        POOLED_SEQTABLE.out.filtered_seqtable
    )

    RENAME_ASVS(
        REMOVE_CHIMERAS.out.nonchim_seqtable,
        POOLED_SEQTABLE.out.filtered_seqtable
    )    

    emit:
    nonchimeric_asvs = RENAME_ASVS.out.nonchimeric_asvs
    seqtable_renamed = RENAME_ASVS.out.seqtable_renamed
    readmap = RENAME_ASVS.out.readmap
    inferred = DADA_INFER.out.inferred
    merged_seqs = POOLED_SEQTABLE.out.merged_seqs
    filtered_seqtable = POOLED_SEQTABLE.out.filtered_seqtable
}

