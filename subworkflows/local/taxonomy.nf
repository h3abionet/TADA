include { DADA2_ASSIGN_TAXA_SPECIES    } from '../../modules/local/assigntaxaspecies'
include { DADA2_TAXTABLE2TEXT          } from '../../modules/local/taxtable2txt'
include { DADA2_TAXFILTER              } from '../../modules/local/dada2taxfilter'

workflow TAXONOMY {

    take:
    readmap      // TODO: this should be modified to simple FASTA
    ref_file     // General use
    species_file // DADA2-specific
    seqtab       // RDS (TODO: should make this CSV or TSV)

    main:
    ch_versions = Channel.empty()
    ch_taxtab_rds = Channel.empty()
    ch_taxmetrics_rds =  Channel.empty()
    ch_readmap_rds = readmap
    ch_seqtab_rds =  seqtab

    DADA2_ASSIGN_TAXA_SPECIES(
        ch_readmap_rds,
        ref_file,
        species_file
    )

    ch_taxtab = DADA2_ASSIGN_TAXA_SPECIES.out.taxtab_rds
    ch_metrics =  DADA2_ASSIGN_TAXA_SPECIES.out.metrics_rds

    if (params.tax_filter) {
        DADA2_TAXFILTER(
            ch_readmap_rds,
            ch_seqtab_rds,
            ch_taxtab,
            ch_metrics
        )
        ch_readmap_rds = DADA2_TAXFILTER.out.readmap_tax_filtered_rds
        ch_seqtab_rds = DADA2_TAXFILTER.out.seqtab_tax_filtered_rds
        ch_taxtab_rds = DADA2_TAXFILTER.out.taxtab_tax_filtered_rds
        ch_taxmetrics_rds = DADA2_TAXFILTER.out.taxmetrics_tax_filtered_rds
    } 

    emit:
    ch_taxtab_rds = ch_taxtab_rds
    ch_taxmetrics_rds = ch_taxmetrics_rds
    ch_readmap_rds = ch_readmap_rds
    ch_seqtab_rds = ch_seqtab_rds
    versions = ch_versions                        // channel: [ versions.yml ]
}
