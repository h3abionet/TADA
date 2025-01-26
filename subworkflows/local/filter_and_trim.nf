include { ILLUMINA_DADA2_FILTER_AND_TRIM   } from '../../modules/local/illumina_filterandtrim'
include { PACBIO_DADA2_FILTER_AND_TRIM     } from '../../modules/local/pacbio_filterandtrim'
include { PACBIO_CUTADAPT                  } from '../../modules/local/pacbio_cutadapt'
include { MERGE_TRIM_TABLES                } from '../../modules/local/mergetrimtables'

workflow FILTER_AND_TRIM {

    take:
    input           //channel: [val(meta), path(reads)]
    skip_filtering  

    main:
    // Three options for Illumina data:
    //       DADA2 trimming and filtering (PE and SE currently) - implemented
    //       cutadapt-based (primers + Ns) + vsearch (EE) - NYI
    //       Hybrid (variable length) - NYI
    // Two options for PacBio:
    //       DADA2 filtering (filter) - NYI
    //       cutadapt (trim) + vsearch - NYI

    ch_reports = Channel.empty()
    ch_trimmed = Channel.empty()
    ch_trimmed_R1 = Channel.empty()
    ch_trimmed_R2 = Channel.empty()

    // TODO: we're probably going to move to requiring the primer sequences to
    //       make the workflow more flexible re: trimming options, esp. since
    //       the current version assumes the presence of primer sequences and
    //       does a hard trim. This also allows for passing in cutadapt anchors 
    //       and primer options (would need to parse these out)
    for_primer = params.for_primer
    for_primer_rc = ""
    rev_primer = params.rev_primer
    rev_primer_rc = ""

    if (for_primer && rev_primer) {
        for_primer_rc = reverse_complement(for_primer)
        rev_primer_rc = reverse_complement(rev_primer)
    }

    if (params.platform == "pacbio") {

        // TODO: this could be modified/split into a `cutadapt`-only step; there
        // are additional filters for max EE and max N in cutadapt
        PACBIO_CUTADAPT(
            input
        )

        // TODO: should be summarized as well, go to MultiQC
        // ch_reports = PACBIO_CUTADAPT_FILTER_AND_TRIM.out.cutadapt_report.collect()
        
        // TODO: this could be modified/split into a `DADA2`-only step
        PACBIO_DADA2_FILTER_AND_TRIM(
            PACBIO_CUTADAPT.out.cutadapt_trimmed
        )
        ch_trimmed = PACBIO_DADA2_FILTER_AND_TRIM.out.trimmed
        ch_trimmed_R1 = PACBIO_DADA2_FILTER_AND_TRIM.out.trimmed
        ch_reports = PACBIO_DADA2_FILTER_AND_TRIM.out.trimmed_report.collect()
    } else {
        // this handles both paired and single-end data
        ILLUMINA_DADA2_FILTER_AND_TRIM(
            input
        )
        ch_trimmed = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed
        ch_reports = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_report.collect()
        ch_trimmed_R1 = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_R1
        ch_trimmed_R2 = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_R2
    }

    // TODO: add variable-length and PacBio
    MERGE_TRIM_TABLES(
        ch_reports
    )

    // Channel setup

    // We need to group data depending on which downstream steps are needed.  There
    // are two combinations possible

    // 1. The immediate downstream QC steps can use the meta info and the read pairs.
    //    Instead of doing handstands reusing the two channels above, we emit channels 
    //    with the reads paired if needed.

    // 2. LearnErrors and the pooled denoising branch requires all R1 and all R2, but 
    //    the two groups can be processed in parallel.  So we set up the channels with 
    //    this in mind. No sample ID info is really needed.
    ch_trimmed_infer = ch_trimmed_R1
            .map { [ 'R1', it[1]] }
            .concat(ch_trimmed_R2.map {['R2', it[1]] } )
            .groupTuple(sort: true)

    emit:
    trimmed = ch_trimmed
    trimmed_report = MERGE_TRIM_TABLES.out.trimmed_report // channel: [ RDS ]
    trimmed_infer = ch_trimmed_infer
}

// def clean_primers(primer) {
//     // returns a clean primer string, IUPAC codes 
//     // w/o any metadata or anchors. Assumes cutadapt 
//     // filtering
// }

def reverse_complement(primer) {
    // returns the revcomp, handles IUPAC ambig codes
    // tr "[ATGCUNYRSWKMBDHV]" "[TACGANRYSWMKVHDB]"
    return primer.reverse().collect { 
        switch (it) {
            case 'A': return 'T'
            case 'T': return 'A'
            case 'G': return 'C'
            case 'C': return 'G'
            case 'U': return 'A'
            case 'N': return 'N'
            case 'Y': return 'R'
            case 'R': return 'Y'
            case 'S': return 'S'
            case 'W': return 'W'
            case 'K': return 'M'
            case 'M': return 'K'
            case 'B': return 'V'
            case 'D': return 'H'
            case 'H': return 'D'
            case 'V': return 'B'
            default: return it // handle invalid characters if needed
        }
    }.join('')
}
