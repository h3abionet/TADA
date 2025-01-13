include { BIOM                   } from '../../modules/local/biom'
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
    seq_table_rds
    seq_table_qiime
    tax_table_rds
    tax_table_tsv
    asvs
    alignment
    unrooted_tree
    rooted_tree

    main:
    ch_versions = Channel.empty()

    if (params.to_BIOM) {
        BIOM(
            seq_table_rds,
            tax_table_rds
        )
    }

    if (params.to_QIIME2) {
        QIIME2_FEATURETABLE(
            seq_table_qiime
        )

        QIIME2_TAXTABLE(
            tax_table_tsv
        )

        QIIME2_SEQUENCE(
            asvs
        )

        if (!params.skip_alignment) {
            QIIME2_ALIGNMENT(
                alignment
            )
        }

        if (!params.skip_tree) {
            QIIME2_TREE(
                unrooted_tree,
                rooted_tree
            )
        }
    }

    // TODO: May become redundant with versions
    SESSION_INFO()

    emit:
    versions = ch_versions                     // channel: [ versions.yml ]
}

