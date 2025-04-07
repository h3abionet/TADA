#!/usr/bin/env Rscript --vanilla --slave
packages <- c('dada2', 'optparse', 'yaml', 'stringr')

suppressPackageStartupMessages({
    lapply(packages, require, character.only = TRUE)
    })

coerce_to_numeric <- function(option, flag, option_value, parser) {
    print(option)
    # 42
    # data <- str_split_1(option_value, ",")
}

coerce_to_bool <- function(x) {
    
}

# https://stackoverflow.com/questions/9129673/passing-list-of-named-parameters-to-function
option_list = list(
    # specific to this workflow
    make_option("--fwd",
        type="character",
        help="The file path to the fastq file. Compressed file formats such as 
        .fastq.gz and .fastq.bz2 are supported"),
    make_option("--filt", 
        type="character",
        help="The path to the output filtered file corresponding to the fwd 
        input file"),
    make_option("--rev",
        type="character",
        help="The file path to the reverse fastq file from paired-end 
        sequence data corresponding to those provided to the fwd argument"),
    make_option("--filt.rev", 
        type="character",
        help="The path to the output filtered file corresponding to the rev 
        input file"),
    make_option("--multithread", 
        help="If TRUE, input files are filtered in parallel via mclapply. 
        If an integer is provided, it is passed to the mc.cores argument 
        of mclapply. Note that the parallelization here is by forking, 
        and each process is loading another fastq file into memory."),
    make_option("--trimLeft", 
        action="callback",
        callback=coerce_to_numeric,
        default="20,20",
        help="The number of nucleotides to remove from the start of each read. 
        If both truncLen and trimLeft are provided, filtered reads will have 
        length truncLen-trimLeft 
        (use comma-separated list for paired reads)"),
    make_option("--trimRight", 
        help="The number of nucleotides to remove from the end of each read. 
        If both truncLen and trimRight are provided, truncation will be 
        performed after trimRight is enforced 
        (use comma-separated list for paired reads)"),
    make_option(c("--truncLen"), 
        help="Truncate reads after truncLen bases. Reads shorter than this 
        are discarded 
        (use comma-separated list for paired reads)"),
    make_option(c("--truncQ"), 
        help="Truncate reads at the first instance of a quality score less 
        than or equal to this 
        (use comma-separated list for paired reads)"),
    make_option(c("--minQ"),
        help="After truncation, reads contain a quality score less than minQ 
        will be discarded.
        (use comma-separated list for paired reads)"),
    make_option(c("--maxEE"), 
        help="After truncation, reads with higher than maxEE 'expected 
        errors' will be discarded. Expected errors are calculated from the 
        nominal definition of the quality score: 
        EE = sum(10^(-Q/10)) 
        (use comma-separated list for paired reads)"),
    make_option(c("--maxN"), 
        help="Maximum Ns allowed 
        (use comma-separated list to differentiate R1/R2 with paired reads)"),
    make_option(c("--maxLen"), 
        help="Remove reads with length greater than maxLen. maxLen is enforced 
        before trimming and truncation
        (use comma-separated list to differentiate R1/R2 with paired reads)"),
    make_option(c("--minLen"), 
        help="Remove reads with length less than minLen. minLen is enforced 
        after trimming and truncation.
        (use comma-separated list to differentiate R1/R2 with paired reads)"),
    make_option(c("--rm.phix"), 
        help="Discard reads that match against the phiX genome 
        (use comma-separated list to differentiate R1/R2 with paired reads)"),
    make_option(c("--rm.lowcomplex"), 
        help="Reads with an effective number of kmers less than this value 
        will be removed. The effective number of kmers is determined by 
        seqComplexity using a Shannon information approximation. The default
        kmer-size is 2, and therefore perfectly random sequences will 
        approach an effective kmer number of:
        16 = 4 (nucleotides) ^ 2 (kmer size) 
        (use comma-separated list to differentiate R1/R2 with paired reads)"),
    make_option(c("--orient.fwd"), 
        help="A character string present at the start of valid reads. 
        Only allows unambiguous nucleotides. This string is compared to the 
        start of each read, and the reverse complement of each read. If 
        it exactly matches the start of the read, the read is kept. If it 
        exactly matches the start of the reverse-complement read, the read 
        is reverse-complemented and kept. Otherwise the read if filtered 
        out."),
    make_option(c("--matchIDs"), 
        help="Paired-read filtering only. Whether to enforce matching between 
        the id-line sequence identifiers of the forward and reverse 
        fastq files."),
    make_option(c("--n"), 
        help="The number of records (reads) to read in and filter at any 
        one time. This controls the peak memory requirement so that very 
        large fastq files are supported.")
    # make_option(c("--trimFor"), type="numeric", , help="trim from 5' for R1"),
    # make_option(c("--trimRev"), type="numeric", , help="trim from 5' for R2"),
    # make_option(c("--truncFor"), type="numeric", , help="truncate R1 to this"),
    # make_option(c("--truncRev"), type="numeric", , help="trucate R2 to this"),
    # make_option(c("--maxEEFor"), type="numeric", , help="maxEE for R1"),
    # make_option(c("--maxEERev"), type="numeric", , help="maxEE for R2"),
)

required <- c("id", "fwd")

opt <- parse_args(OptionParser(option_list=option_list, formatter=TitledHelpFormatter))

# print(opt)

# tmpFunc <- function(...) {
#     print(as.list(match.call()))
#     filterAndTrim(...)
# }

# TODO switch to a do.call() using input args
# out <- filterAndTrim(fwd        = opt$fwd,
#                     filt        = paste0(opt$id, ".R1.filtered.fastq.gz"),
#                     rev         = opt$rev,
#                     filt.rev    = paste0(opt$id, ".R2.filtered.fastq.gz"),
#                     trimLeft    = c(opt$trimFor,  opt$trimRev),  # fix
#                     truncLen    = c(opt$truncFor, opt$truncRev), # fix
#                     maxEE       = c(opt$maxEEFor, opt$maxEERev), # fix
#                     truncQ      = opt$truncQ,
#                     maxN        = opt$maxN,
#                     rm.phix     = opt$rm.phix,
#                     maxLen      = opt$maxLen,
#                     minLen      = opt$minLen,
#                     compress    = TRUE,
#                     verbose     = TRUE,
#                     multithread = opt$cpus
#                     )

# colnames(out) <- c('input', 'filtered')

# write.csv(out, paste0(opt$id, ".trimmed.txt"))

# out <- filterAndTrim(fwd        = "${reads[0]}",
#                     filt        = "${meta.id}.R1.filtered.fastq.gz",
#                     rev         = if("${reads[1]}" == "null") NULL else "${reads[1]}",
#                     filt.rev    = if("${reads[1]}" == "null") NULL else "${meta.id}.R2.filtered.fastq.gz",
#                     trimLeft    = if("${reads[1]}" == "null") ${params.trim_for} else  c(${params.trim_for}, ${params.trim_rev}),
#                     truncLen    = if("${reads[1]}" == "null") ${params.trunc_for} else c(${params.trunc_for}, ${params.trunc_rev}),
#                     maxEE       = if("${reads[1]}" == "null") ${params.maxEE_for} else c(${params.maxEE_for}, ${params.maxEE_rev}), 
#                     truncQ      = ${params.truncQ},
#                     maxN        = ${params.maxN},
#                     rm.phix     = as.logical("${params.rmPhiX}"),
#                     maxLen      = ${params.max_read_len},
#                     minLen      = ${params.min_read_len},
#                     compress    = TRUE,
#                     verbose     = TRUE,
#                     multithread = ${task.cpus}
#                     )

# colnames(out) <- c('input', 'filtered')

# write.csv(out, "${meta.id}.trimmed.txt")

# if (opt$session) {
#     pv <- sapply(packages, function(x) { as.character(packageVersion(x)) }, simplify = FALSE)
#     pv$R <- paste0(R.Version()[c("major", "minor")], collapse = ".")
#     session_list <- list(pv)
#     names(session_list) <- opt$process
#     write_yaml(session_list, "versions.yml")
# }
