// include { MMSEQS_CONVERTALIS } from '../../../modules/local/mmseqs/convertalis'
// include { MMSEQS_CONVERTMSA  } from '../../../modules/local/mmseqs/convertmsa'
// include { MMSEQS_MSA2PROFILE } from '../../../modules/local/mmseqs/msa2profile'
// include { MMSEQS_SEARCH      } from '../../../modules/nf-core/mmseqs/search'
// include { MMSEQS_CREATEDB    as MMSEQS_CREATEDB_PROFILE } from '../../../modules/nf-core/mmseqs/createdb'
// include { MMSEQS_CREATEDB    as MMSEQS_CREATEDB_QUERY   } from '../../../modules/nf-core/mmseqs/createdb'
// include { MMSEQS_FILTER_DATA }

// TODO: split up and move into modules 
process MMSEQS_PROFILE_FULL {
    tag "${profile.getSimpleName()}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    // tuple val(meta), path("${prefix}/"), emit: db_search
    path(profile)
    path(asvs)
    path(readmap)

    output:
    path("asv_vs_${profile.getSimpleName()}_filtering_results.profile.m8"), emit: db_search
    path("asv_vs_${profile.getSimpleName()}_filtering_results.profile.ids"), emit: db_ids
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    ## TODO: split into separate steps!!!
    mmseqs convertmsa ${profile} \\
        ${profile.getSimpleName()}.msa_db
    
    mmseqs msa2profile \\
        ${profile.getSimpleName()}.msa_db \\
        ${profile.getSimpleName()}.profile \\
        --match-mode 1 --threads ${task.cpus}

    mmseqs createindex \\
        ${profile.getSimpleName()}.profile \\
        tmp \\
        --threads ${task.cpus} \\        
        -k 6 -s 7

    mmseqs createdb \\
        ${asvs} \\
        asvs \\
        --dbtype 2

    mkdir ${profile.getSimpleName()}

    # most sensitive setting for matches
    mmseqs search \\
        asvs \\
        ${profile.getSimpleName()}.profile \\
        ${profile.getSimpleName()}/asvs.results  \\
        tmp \\
        --threads ${task.cpus} \\
        -k 6 -s 7 

    mmseqs convertalis \\
        asvs \\
        ${profile.getSimpleName()}.profile \\
        asvs.results asv_vs_${profile.getSimpleName()}_filtering_results.profile.m8

    cut -f1 asv_vs_${profile.getSimpleName()}_filtering_results.profile.m8 | \\
        sort | \\
        uniq > \\
        asv_vs_${profile.getSimpleName()}_filtering_results.profile.ids

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}

process MMSEQS_DATABASE_FULL {
    tag "${database.getSimpleName()}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2:17.b804f--hd6d6fdc_1':
        'biocontainers/mmseqs2:17.b804f--hd6d6fdc_1' }"

    input:
    // tuple val(meta), path("${prefix}/"), emit: db_search
    path(database)
    path(asvs)

    output:
    path("${asvs.getSimpleName()}_vs_${database.getSimpleName()}.database.m8"), emit: db_search
    path("${asvs.getSimpleName()}_vs_${database.getSimpleName()}.database.ids"), emit: db_ids
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    # most sensitive setting for matches
    mmseqs easy-search ${asvs} \\
        ${database} \\
        ${asvs.getSimpleName()}_vs_${database.getSimpleName()}.database.m8 \\
        tmp \\
        -s 7.5 \\
        --threads ${task.cpus}
    
    cut -f1 ${asvs.getSimpleName()}_vs_${database.getSimpleName()}.database.m8 | \\
        sort | uniq > \\
        ${asvs.getSimpleName()}_vs_${database.getSimpleName()}.database.ids

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(mmseqs | grep 'Version' | sed 's/MMseqs2 Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    touch ${prefix}.

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mmseqs: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
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
    path("asvs.${params.id_type}.nochim.filtered.fna"), emit: filtered_asvs
    path("seqtab_final.${params.id_type}.filtered.RDS"), emit: filtered_seqtab
    path("readmap.${params.id_type}.filtered.RDS"), emit: filtered_readmap

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
    seqtab <- readRDS("${seqtab}")
    readmap <- readRDS("${readmap}")
    asvs <- readDNAStringSet("${asvs}", format="fasta")
    ids <- readLines("${ids}")
    writeXStringSet(asvs[ids], "asvs.${params.id_type}.nochim.filtered.fna")
    saveRDS(seqtab[,ids], "seqtab_final.${params.id_type}.filtered.RDS")
    saveRDS(readmap[readmap\$ids %in% ids,], "readmap.${params.id_type}.filtered.RDS")
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

    if (params.filter == "mmseqs_aa_profile") {

        // TODO: these need to be sanity checked early on
        ch_profile = file(params.mmseqs_profile, checkIfExists: true)

        // TODO: we need to allow more parameters through
        MMSEQS_PROFILE_FULL(
            ch_profile,
            ch_asvs
        )
        ch_versions = ch_versions.mix(MMSEQS_PROFILE_FULL.out.versions)
        ch_filtered_ids = MMSEQS_PROFILE_FULL.out.db_ids
    } else if (params.filter == "mmseqs_aa_database") {

        // TODO: these need to be sanity checked early on
        ch_database = file(params.mmseqs_database, checkIfExists: true)

        // TODO: we need to allow more parameters through   
        MMSEQS_DATABASE_FULL(
            ch_database,
            ch_asvs
        )

        ch_versions = ch_versions.mix(MMSEQS_DATABASE_FULL.out.versions)
        ch_filtered_ids = MMSEQS_DATABASE_FULL.out.db_ids
    }

    FILTER_TADA_DATA(
        ch_seqtab,
        ch_asvs,
        ch_readmap,
        ch_filtered_ids
    )

    emit:
    ch_filtered_seqtab = params.filter_dryrun ? 
        ch_seqtab :
        FILTER_TADA_DATA.out.filtered_seqtab
    ch_filtered_asvs = params.filter_dryrun ? 
        ch_asvs :
        FILTER_TADA_DATA.out.filtered_asvs
    ch_filtered_readmap = params.filter_dryrun ? 
        ch_readmap :
        FILTER_TADA_DATA.out.filtered_readmap
    versions = ch_versions // channel: [ versions.yml ]
}

