include { DADA2_ASSIGN_TAXA_SPECIES    } from '../../modules/local/assigntaxaspecies'
include { DADA2_TAXTABLE2TEXT          } from '../../modules/local/taxtable2txt'
// include { DADA2_TAXFILTER              } from '../../modules/local/dada2taxfilter'

workflow TAXONOMY {

    take:
    readmap         // TODO: this should be modified to simple FASTA
    ref_file     // General use
    species_file // DADA2-specific
    // seqtab       // RDS (TODO: should make this CSV or TSV)

    main:
    ch_versions = Channel.empty()
    ch_taxtab = Channel.empty()
    ch_metrics =  Channel.empty()
    // ch_filtered_seqtab =  Channel.empty()
    // TODO: eventually this will have multiple options for
    //       taxonomic assignment
    DADA2_ASSIGN_TAXA_SPECIES(
        readmap,
        ref_file,
        species_file
    )

    // harmonize to TSV tax assn table + any metrics per ASV ID
    DADA2_TAXTABLE2TEXT(
        DADA2_ASSIGN_TAXA_SPECIES.out.taxtab_rds,
        DADA2_ASSIGN_TAXA_SPECIES.out.metrics_rds
    )

    // if (params.tax_filter) {
    //     DADA2_TAXFILTER(
    //         readmap,
    //         DADA2_ASSIGN_TAXA_SPECIES.out.taxtab_rds,
    //         asvs,
    //         seqtab
    //     )
    //     ch_taxtab = 
    //     ch_metrics =  Channel.empty()
    //     ch_filtered_seqtab = seqtab
    // } else {

    // }

    emit:
    // TODO: I'd like to make the data less dependent on R, so 
    // move to passing standard formats (TSV, CSV, FASTA, etc)
    ch_taxtab_rds = DADA2_ASSIGN_TAXA_SPECIES.out.taxtab_rds 
    ch_taxtab = DADA2_TAXTABLE2TEXT.out.taxtab    // channel: [ TSV ]
    ch_metrics = DADA2_TAXTABLE2TEXT.out.metrics  // channel: [ TSV ]
    // ch_seqtab = DADA2_TAXTABLE2TEXT.out.seqtab  // channel: [ TSV ]
    versions = ch_versions                        // channel: [ versions.yml ]
}
