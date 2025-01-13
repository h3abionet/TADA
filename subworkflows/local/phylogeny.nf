include { DECIPHER               } from '../../modules/local/decipher'
include { PHANGORN               } from '../../modules/local/phangorn'
include { FASTTREE               } from '../../modules/local/fasttree'
include { ROOT_TREE              } from '../../modules/local/roottree'

workflow PHYLOGENY {
    take:
    asvs

    main:

    ch_alignment = Channel.empty()
    ch_unrooted_tree = Channel.empty()
    ch_rooted_tree = Channel.empty()
    ch_versions = Channel.empty()

    if (!params.skip_alignment) {
        DECIPHER(
            asvs
        )
        ch_alignment = DECIPHER.out.alignment
        if (!params.skip_tree) {
            if (params.phylo_tool == 'phangorn') {
                PHANGORN(
                    ch_alignment
                )
                ch_unrooted_tree = PHANGORN.out.treeGTR
            } else if (params.phylo_tool == 'fasttree') {
                FASTTREE(
                    ch_alignment
                )
                ch_unrooted_tree = FASTTREE.out.treeGTR
            }

            ROOT_TREE(
                ch_unrooted_tree,
                params.phylo_tool
            )
            ch_rooted_tree = ROOT_TREE.out.rooted_tree
        }
    }

    emit:
    ch_alignment
    ch_unrooted_tree
    ch_rooted_tree
    versions = ch_versions // channel: [ versions.yml ]
}

