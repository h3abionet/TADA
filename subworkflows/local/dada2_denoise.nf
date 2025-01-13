// TODO: move to a subworkflow and implement pooled vs per-sample + optional priors
include { ILLUMINA_DADA2_LEARN_ERRORS           } from '../../modules/local/illumina_learnerrors'
include { PACBIO_DADA2_LEARN_ERRORS             } from '../../modules/local/pacbio_learnerrors'
include { DADA_INFER                            } from '../../modules/local/dadainfer'
include { POOLED_SEQTABLE                       } from '../../modules/local/pooledseqtable'
include { DADA2_REMOVE_CHIMERAS                 } from '../../modules/local/removechimeras'
include { RENAME_ASVS                           } from '../../modules/local/renameasvs'
include { DADA2_SEQTABLE2TEXT                   } from '../../modules/local/seqtable2txt'

workflow DADA2_DENOISE {

    take:
    ch_trimmed_infer 

    main:

    ch_versions = Channel.empty()

    ch_infer = Channel.empty()
    // TODO nf-core: substitute modules here for the modules of your subworkflow
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

    DADA2_REMOVE_CHIMERAS(
        POOLED_SEQTABLE.out.filtered_seqtable
    )

    RENAME_ASVS(
        DADA2_REMOVE_CHIMERAS.out.nonchim_seqtable,
        POOLED_SEQTABLE.out.filtered_seqtable
    )

    DADA2_SEQTABLE2TEXT(
        RENAME_ASVS.out.seqtable_renamed
    )

    emit:
    seqtab2qiime = DADA2_SEQTABLE2TEXT.out.seqtab2qiime
    nonchimeric_asvs = RENAME_ASVS.out.nonchimeric_asvs
    seqtable_renamed = RENAME_ASVS.out.seqtable_renamed
    readmap = RENAME_ASVS.out.readmap
    inferred = DADA_INFER.out.inferred
    merged_seqs = POOLED_SEQTABLE.out.merged_seqs
    filtered_seqtable = POOLED_SEQTABLE.out.filtered_seqtable
}

