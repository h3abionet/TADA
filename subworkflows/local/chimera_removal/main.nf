include { DADA2_REMOVE_CHIMERAS                 } from '../../../modules/local/removechimeras'
include { RENAME_ASVS                           } from '../../../modules/local/renameasvs'

workflow CHIMERA_REMOVAL {

    take:
    filtered_seqtable // channel: RDS

    main:

    ch_versions = Channel.empty()

    DADA2_REMOVE_CHIMERAS(
        filtered_seqtable
    )

    RENAME_ASVS(
        DADA2_REMOVE_CHIMERAS.out.nonchim_seqtable,
        filtered_seqtable
    )

    emit:
    nonchimeric_asvs = RENAME_ASVS.out.nonchimeric_asvs
    seqtable_renamed = RENAME_ASVS.out.seqtable_renamed
    readmap = RENAME_ASVS.out.readmap
    versions = ch_versions                     // channel: [ versions.yml ]
}

