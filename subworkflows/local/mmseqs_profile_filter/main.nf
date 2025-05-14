include { MMSEQS_CREATEDB       } from '../../../modules/nf-core/mmseqs/createdb'
include { MMSEQS_EASYSEARCH     } from '../../../modules/nf-core/mmseqs/easysearch/main'

// process MMSEQS_PROFILE_FULL {
//     tag "${profile.getSimpleName()}"
//     label 'process_low'

//     container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
//         'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
//         'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

//     input:
//     // tuple val(meta), path("${prefix}/"), emit: db_search
//     path(profile)
//     path(asvs)
//     path(readmap)

//     output:
//     path("asv_vs_${profile.getSimpleName()}_filtering_results.profile.m8"), emit: db_search
//     path("asv_vs_${profile.getSimpleName()}_filtering_results.profile.ids"), emit: db_ids
//     path "versions.yml"           , emit: versions

//     when:
//     task.ext.when == null || task.ext.when

//     script:
//     def args = task.ext.args ?: ''
    
//     """
//     ## TODO: split into separate steps!!!
//     mmseqs convertmsa ${profile} \\
//         ${profile.getSimpleName()}.msa_db
    
//     mmseqs msa2profile \\
//         ${profile.getSimpleName()}.msa_db \\
//         ${profile.getSimpleName()}.profile \\
//         --match-mode 1 --threads ${task.cpus}

//     mmseqs createindex \\
//         ${profile.getSimpleName()}.profile \\
//         tmp \\
//         --threads ${task.cpus} \\        
//         -k 6 -s 7

//     mmseqs createdb \\
//         ${asvs} \\
//         asvs \\
//         --dbtype 2

//     mkdir ${profile.getSimpleName()}

//     # most sensitive setting for matches
//     mmseqs search \\
//         asvs \\
//         ${profile.getSimpleName()}.profile \\
//         ${profile.getSimpleName()}/asvs.results  \\
//         tmp \\
//         --threads ${task.cpus} \\
//         -k 6 -s 7 

//     mmseqs convertalis \\
//         asvs \\
//         ${profile.getSimpleName()}.profile \\
//         asvs.results asv_vs_${profile.getSimpleName()}_filtering_results.profile.m8

//     cut -f1 asv_vs_${profile.getSimpleName()}_filtering_results.profile.m8 | \\
//         sort | \\
//         uniq > \\
//         asv_vs_${profile.getSimpleName()}_filtering_results.profile.ids

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
//     END_VERSIONS
//     """

//     stub:
//     def args = task.ext.args ?: ''
    
//     """
//     touch ${prefix}.bam

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         mmseqs: \$(samtools --version |& sed '1!d ; s/samtools //')
//     END_VERSIONS
//     """
// }

// TODO: this can be merged in the R-based FILTER_TADA_DATA step
process MMSEQS_DATABASE_FILTER {
    tag "${meta2.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    tuple val(meta), path(tsv)
    tuple val(meta2), path(db)

    output:
    path("asvs_vs_${meta2.id}.database.ids"), emit: db_ids
    // path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """    
    cut -f1 ${tsv} | \\
        sort | uniq > \\
        asvs_vs_${meta2.id}.database.ids
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch asvs_vs_${prefix}.database.ids
    """
}

process FILTER_TADA_DATA {
    tag "FILTER_TADA_DATA"
    label 'process_low'

    container "ghcr.io/h3abionet/tada:dev"

    input:
    path(seqtab)
    path(asvs)
    path(readmap)
    file(ids)

    output:
    path("asvs.${params.id_type}.search_filtered.fna"), emit: filtered_asvs
    path("seqtab.${params.id_type}.search_filtered.RDS"), emit: filtered_seqtab
    path("readmap.${params.id_type}.search_filtered.RDS"), emit: filtered_readmap

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    """
    #!/usr/bin/env Rscript

    # TODO: clean me up!
    suppressPackageStartupMessages(library(dada2))
    suppressPackageStartupMessages(library(Biostrings))
    suppressPackageStartupMessages(library(tidyverse))
    #mmseqs_results <- read_tsv("", 
    #    col_names=c(
    #    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    #    "qstart", "qend", "tstart", "tend", "evalue", "bits"))

    seqtab <- readRDS("${seqtab}")
    readmap <- readRDS("${readmap}")
    asvs <- readDNAStringSet("${asvs}", format="fasta")
    ids <- readLines("${ids}")
    writeXStringSet(asvs[ids], "asvs.${params.id_type}.search_filtered.fna")
    saveRDS(seqtab[,ids, drop=FALSE], "seqtab.${params.id_type}.search_filtered.RDS")
    saveRDS(readmap[readmap\$id %in% ids, ,drop=FALSE], "readmap.${params.id_type}.search_filtered.RDS")
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix
    """
    # TODO: make a proper stub
    """
}

workflow MMSEQS_FILTER {

    take:
    ch_seqtab  
    ch_asvs
    ch_readmap
    
    main:
    ch_versions = Channel.empty()
    ch_filtered_ids = Channel.empty()

    // there is no meta here, so we create one
    ch_asvs_meta = ch_asvs
        .map { it -> tuple([id: "asvs"], it)}

    // TODO: sanity check this
    if (params.mmseqs_method == "search") {

        ch_database = Channel.empty()
        
        if (params.mmseqs_fasta) {
            // Create a meta file on the fly
            ch_mmseqs_fasta = Channel
                .fromPath(params.mmseqs_fasta, checkIfExists: true)
                .map { it -> tuple( [id: it.getSimpleName()], it) }
                .first()

            MMSEQS_CREATEDB(
                ch_mmseqs_fasta
            )

            ch_database = MMSEQS_CREATEDB.out.db
            ch_versions = ch_versions.mix(MMSEQS_CREATEDB.out.versions)
        } else {
            ch_database = Channel
                .fromPath(params.mmseqs_database, checkIfExists: true)
                .map { it -> tuple( [id: it.getSimpleName()], it) }
                .first()
        }       

        MMSEQS_EASYSEARCH(
            ch_asvs_meta,
            ch_database
        )
        ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions)

        // TODO: we need to allow more parameters through   
        MMSEQS_DATABASE_FILTER(
            MMSEQS_EASYSEARCH.out.tsv,
            ch_database
        )

        ch_filtered_ids = MMSEQS_DATABASE_FILTER.out.db_ids
    } 

    FILTER_TADA_DATA(
        ch_seqtab,
        ch_asvs,
        ch_readmap,
        ch_filtered_ids
    )

    emit:
    ch_filtered_seqtab = params.search_filter_dryrun ? 
        ch_seqtab :
        FILTER_TADA_DATA.out.filtered_seqtab
    ch_filtered_asvs = params.search_filter_dryrun ? 
        ch_asvs :
        FILTER_TADA_DATA.out.filtered_asvs
    ch_filtered_readmap = params.search_filter_dryrun ? 
        ch_readmap :
        FILTER_TADA_DATA.out.filtered_readmap
    versions = ch_versions // channel: [ versions.yml ]
}

