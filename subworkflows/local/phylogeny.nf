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
        ch_alignment = Channel.empty()
        // This is needed b/c nf-core modules have a 
        // meta file, so we need to create a simple 
        // fake one for the workflow
        // ch_aln = asvs
        //     .map{
        //         [id: params.aligner ], it ->
        //             [ meta, fastafile ]
        //     }
        // if (params.aligner == "decipher") {
            DECIPHER(
                asvs
            )
            ch_alignment = DECIPHER.out.alignment
        // } else if (params.aligner == "mafft") {
        //     MAFFT_ALIGN (
        //         ch_aln ,
        //         [ [:], [] ],
        //         [ [:], [] ],
        //         [ [:], [] ],
        //         [ [:], [] ],
        //         [ [:], [] ],
        //         true
        //     )
        //     ch_alignment = MAFFT_ALIGN.out.fas // the MAFFT module calls its output fas instead of alignment
        //     ch_versions = ch_versions.mix(MAFFT_ALIGN.out.versions.first())
        // } else if (params.aligner == "muscle") {
        //     MUSCLE5_ALIGN (
        //         ch_aln ,
        //         true
        //     )
        //     ch_alignment = MUSCLE5_ALIGN.out.alignment.first() // the MAFFT module calls its output fas instead of alignment
        //     ch_versions = ch_versions.mix(MUSCLE5_ALIGN.out.versions.first())
        // }
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

