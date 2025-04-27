include { BIOM                   } from '../../modules/local/biom'
include { DADA2_SEQTABLE2TEXT    } from '../../modules/local/seqtable2txt'
include { DADA2_TAXTABLE2TEXT    } from '../../modules/local/taxtable2txt'
include { QIIME2_FEATURETABLE    } from '../../modules/local/qiime2featuretable'
include { QIIME2_TAXTABLE        } from '../../modules/local/qiime2taxtable'
include { QIIME2_SEQUENCE        } from '../../modules/local/qiime2seqs'
include { QIIME2_ALIGNMENT       } from '../../modules/local/qiime2aln'
include { QIIME2_TREE            } from '../../modules/local/qiime2tree'
include { SESSION_INFO           } from '../../modules/local/rsessioninfo'

workflow GENERATE_OUTPUT {

    // TODO: I'd like to have this simply be TSV files (no RDS)
    //       so we can generate from other subworkflows if needed
    take:
    seqtable_rds
    // seq_table_qiime
    taxtable_rds
    metrics_rds
    asvs
    alignment
    unrooted_tree
    rooted_tree

    main:
    ch_versions = Channel.empty()
    ch_taxtable_tsv = Channel.empty()

    DADA2_SEQTABLE2TEXT(
        seqtable_rds
    )

    if (params.reference) {
        DADA2_TAXTABLE2TEXT(
            taxtable_rds, metrics_rds
        )
        ch_taxtable_tsv = DADA2_TAXTABLE2TEXT.out.taxtab
    }

    if (params.to_BIOM) {
        BIOM(
            seqtable_rds,
            taxtable_rds
        )
    }

    if (params.to_QIIME2) {

        QIIME2_FEATURETABLE(
            DADA2_SEQTABLE2TEXT.out.seqtab2qiime
        )

        ch_versions = ch_versions.mix(QIIME2_FEATURETABLE.out.versions)

        QIIME2_TAXTABLE(
            ch_taxtable_tsv
        )

        ch_versions = ch_versions.mix(QIIME2_TAXTABLE.out.versions)

        QIIME2_SEQUENCE(
            asvs
        )

        ch_versions = ch_versions.mix(QIIME2_SEQUENCE.out.versions)

        if (!params.skip_alignment) {
            QIIME2_ALIGNMENT(
                alignment
            )
            ch_versions = ch_versions.mix(QIIME2_ALIGNMENT.out.versions)
        }

        if (!params.skip_tree) {
            QIIME2_TREE(
                unrooted_tree,
                rooted_tree
            )

            ch_versions = ch_versions.mix(QIIME2_TREE.out.versions)
        }
    }

    // TODO: May become redundant with versions
    SESSION_INFO()

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}

