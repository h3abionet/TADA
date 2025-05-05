include { DADA2_POOLED_INFER                    } from '../../../modules/local/dadainfer'
include { DADA2_POOLED_SEQTABLE                 } from '../../../modules/local/pooledseqtable'

workflow DADA2_POOLED_DENOISE {

    take:
    ch_infer
    ch_trimmed_infer

    main:
    ch_versions = Channel.empty()
    
    DADA2_POOLED_INFER(ch_infer)

    ch_trimmed = ch_trimmed_infer
        .map { it[1] }
        .flatten()
        .collect()

    DADA2_POOLED_SEQTABLE(
        DADA2_POOLED_INFER.out.inferred.collect(),
        ch_trimmed)

    emit:
    inferred = DADA2_POOLED_INFER.out.inferred
    merged_seqs = DADA2_POOLED_SEQTABLE.out.merged_seqs
    filtered_seqtable = DADA2_POOLED_SEQTABLE.out.filtered_seqtable
    versions = ch_versions                     // channel: [ versions.yml ]
}
