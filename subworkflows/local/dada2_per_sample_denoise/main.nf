include { PER_SAMPLE_INFER                         } from '../../../modules/local/persampleinfer'
include { PER_SAMPLE_MERGE                         } from '../../../modules/local/persamplemerge'
include { PER_SAMPLE_SEQTABLE as PER_SAMPLE_SEQTABLE_SE      } from '../../../modules/local/persampleseqtable'
include { PER_SAMPLE_SEQTABLE as PER_SAMPLE_SEQTABLE_PE      } from '../../../modules/local/persampleseqtable'
include { PER_SAMPLE_TRACKING                      } from '../../../modules/local/persampletracking'

workflow DADA2_PER_SAMPLE_DENOISE {

    take:
    ch_dereps    // channel: [ meta, derep.RDS ]
    ch_errs      // channel: [ err.R[1,2].RDS ]
    for_priors   // channel: [ optional: for_priors.fa ]
    rev_priors   // channel: [ optional: rev_priors.fa ]
    stage        // String, simple label for round

    main:
    ch_inferred = Channel.empty()
    ch_merged = Channel.empty()
    ch_filtered_seqtable = Channel.empty()
    ch_for_priors = Channel.empty()
    ch_rev_priors = Channel.empty()
    ch_readtracking = Channel.empty()
    ch_versions = Channel.empty()

    PER_SAMPLE_INFER(
        ch_dereps,
        ch_errs,   
        for_priors, 
        rev_priors, 
        stage
        )

    ch_per_sample_inferred = PER_SAMPLE_INFER.out.dds

    ch_merging = ch_per_sample_inferred
        .filter { meta, dds -> meta.single_end == false } 
        .map{ meta, dds -> [meta.id, meta, dds] }
        .join( ch_dereps.map { meta, derep -> [meta.id, meta, derep] })
        .map{ id, meta1, dds, meta2, derep -> [ meta1, dds, derep ]}

    PER_SAMPLE_MERGE(
        ch_merging,
        stage
    )

    ch_per_sample_merged = PER_SAMPLE_MERGE.out.merged_reads

    // this should technically be renamed
    PER_SAMPLE_TRACKING(
        ch_per_sample_inferred
        .map { it -> it[1] } 
        .collect(),
        stage
    )

    ch_for_priors = PER_SAMPLE_TRACKING.out.for_priors
    ch_rev_priors = PER_SAMPLE_TRACKING.out.rev_priors
    ch_inferred = PER_SAMPLE_TRACKING.out.inferred.collect()
    ch_readtracking = ch_readtracking.mix(PER_SAMPLE_TRACKING.out.readtracking.collect())

    ch_readmode = PER_SAMPLE_INFER.out.readmode.first()

    PER_SAMPLE_SEQTABLE_SE(
        ch_per_sample_inferred
            .filter{ it[0].single_end }
            .map { it[1] }
            .collect(),
        ch_readmode,
        stage
        )

    PER_SAMPLE_SEQTABLE_PE(
        ch_per_sample_merged.map { it[1] }.collect(),
        ch_readmode,
        stage
        )

    // ch_merged = PER_SAMPLE_SEQTABLE_PE.out.merged_seqs
    // TODO: this probably needs a sanity check and error out if the 
    //       channel includes more than one seqtable
    ch_filtered_seqtable = PER_SAMPLE_SEQTABLE_SE.out.filtered_seqtable.mix(PER_SAMPLE_SEQTABLE_PE.out.filtered_seqtable)
    ch_readtracking = ch_readtracking.mix(PER_SAMPLE_SEQTABLE_SE.out.readtracking, PER_SAMPLE_SEQTABLE_PE.out.readtracking)

    emit:
    inferred = ch_inferred                       // channel: [ meta, [ ddsR1.RDS, ddsR2.RDS? ] ]
    merged_seqs = ch_merged                      // channel: [ meta, [ merged.RDS ] ]
    filtered_seqtable = ch_filtered_seqtable     // channel: [ seqtable.RDS ]
    for_priors = ch_for_priors                   // channel: [ priors.fa    ]
    rev_priors = ch_rev_priors                   // channel: [ priors.fa    ]  
    versions = ch_versions                       // channel: [ versions.yml ]
    readtracking = ch_readtracking               // channel: [ data.csv ]
}
