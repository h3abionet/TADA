include { ILLUMINA_DADA2_FILTER_AND_TRIM   } from '../../modules/local/illumina_filterandtrim'
include { PACBIO_DADA2_FILTER_AND_TRIM     } from '../../modules/local/pacbio_filterandtrim'
include { ILLUMINA_CUTADAPT                } from '../../modules/local/illumina_cutadapt'
include { PACBIO_CUTADAPT                  } from '../../modules/local/pacbio_cutadapt'
include { MERGE_TRIM_TABLES                } from '../../modules/local/mergetrimtables'

workflow FILTER_AND_TRIM {

    take:
    input           //channel: [val(meta), path(reads)]
    skip_filtering

    main:
    ch_reports = Channel.empty()
    ch_trimmed = Channel.empty()
    ch_trimmed_R1 = Channel.empty()
    ch_trimmed_R2 = Channel.empty()
    ch_multiqc_files = Channel.empty()

    for_primer_rc = params.for_primer ? reverse_complement(params.for_primer) : ""
    rev_primer_rc = params.rev_primer ? reverse_complement(params.rev_primer) : ""

    // Two modules/subworkflows for Illumina data:
    //       DADA2 trimming and filtering (PE and SE currently) - implemented
    //       cutadapt-based (primers + Ns) + vsearch (EE) - implemented
    // Two modules/subworkflows for PacBio:
    //       DADA2 filtering (filter) - NYI
    //       cutadapt (trim) - implemented

    if (params.platform == "pacbio") {

        PACBIO_CUTADAPT(
            input
            params.for_primer
            rev_primer_rc
        )

        // TODO: this could be modified/split into a `DADA2`-only step
        // PACBIO_DADA2_FILTER_AND_TRIM(
        //     PACBIO_CUTADAPT.out.cutadapt_trimmed
        // )
        ch_trimmed = PACBIO_CUTADAPT.out.trimmed
        ch_trimmed_R1 = PACBIO_CUTADAPT.out.trimmed
        ch_reports = PACBIO_CUTADAPT.out.cutadapt_report.collect()
        ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA_CUTADAPT.out.cutadapt_json)
    } else {
        // this handles both paired and single-end data
        if (params.trimmer == "dada2") {
            ILLUMINA_DADA2_FILTER_AND_TRIM(
                input 
            )
            ch_trimmed = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed
            ch_reports = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_report.collect()
            ch_trimmed_R1 = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_R1
            ch_trimmed_R2 = ILLUMINA_DADA2_FILTER_AND_TRIM.out.trimmed_R2
        } else if (params.trimmer == "cutadapt") {
            ILLUMINA_CUTADAPT(
                input,
                params.for_primer,
                params.rev_primer,
                for_primer_rc,
                rev_primer_rc
            )
            ch_trimmed = ILLUMINA_CUTADAPT.out.trimmed
            ch_reports = ILLUMINA_CUTADAPT.out.trimmed_report.collect()
            ch_trimmed_R1 = ILLUMINA_CUTADAPT.out.trimmed_R1
            ch_trimmed_R2 = ILLUMINA_CUTADAPT.out.trimmed_R2
            ch_multiqc_files = ch_multiqc_files.mix(ILLUMINA_CUTADAPT.out.cutadapt_json)
        }
    }

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
    ch_multiqc_files
}

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

// def clean_primers(primer) {
//     // returns a clean primer string, IUPAC codes 
//     // w/o any metadata or anchors. Assumes cutadapt 
//     // filtering
// }

