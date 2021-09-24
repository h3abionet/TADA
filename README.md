![Github CI Status](https://github.com/grbot/TADA/actions/workflows/ci.yml/badge.svg)
[![Travis-CI Build Status](https://travis-ci.com/h3abionet/TADA.svg?branch=master)](https://travis-ci.com/h3abionet/TADA)
[![DOI](https://zenodo.org/badge/218786496.svg)](https://zenodo.org/badge/latestdoi/218786496)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

<p>
<img align="left" src="./assets/cbio_logo.png" width="240" hspace="50"/>
<img align="left" src="./assets/HPCBio.png" width="350" hspace="50"/>
</br></br></br>
</p>

# TADA - Targeted Amplicon Diversity Analysis using DADA2, implemented in Nextflow

A dada2-based workflow using the Nextflow workflow manager.  The basic pipeline is currently implemented, including some basic read-tracking. This pipeline is adapted from https://github.com/HPCBio/dada2-Nextflow for implementation on the UCT high-performance compute cluster

## Basic usage:

    This pipeline can be run specifying parameters in a config file or with command line flags.

    The typical example for running the pipeline with command line flags is as follows:
    nextflow run h3abionet/TADA --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex

    The typical command for running the pipeline with your own config (instead of command line flags) is as follows:
    nextflow run h3abionet/TADA -c dada2_user_input.config -profile uct_hex
    where:
    dada2_user_input.config is the configuration file (see example 'dada2_user_input.config')
    NB: -profile uct_hex still needs to be specified from the command line

    To override existing values from the command line, please type these parameters:

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Hardware config to use. Currently profile available for UCT's HPC 'uct_hex' - create your own if necessary
      --trimFor                     integer. headcrop of read1 (set 0 if no trimming is needed)
      --trimRev                     integer. headcrop of read2 (set 0 if no trimming is needed)
      --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz)
      --amplicon		                Type of analysis (16S or ITS)
      --runtree                     phangorn or fasttree

    All available read preparation parameters:
      --trimFor                     integer. headcrop of read1
      --trimRev                     integer. headcrop of read2
      --truncFor                    nteger. truncate read1 here (i.e. if you want to trim 10bp off the end of a 250bp R1, truncFor should be set to 240). enforced before trimFor/trimRev
      --truncRev                    nteger. truncate read2 here (i.e. if you want to trim 10bp off the end of a 250bp R2, truncRev should be set to 240). enforced before trimFor/trimRev
      --maxEEFor                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be  discarded. EE = sum(10^(-Q/10)), default=2
      --maxEERev                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be  discarded. EE = sum(10^(-Q/10)), default=2
      --truncQ                      integer. Truncate reads at the first instance of a quality score less than or equal to  truncQ; default=2
      --maxN                        integer. Discard reads with more than maxN number of Ns in read; default=0
      --maxLen                      integer. maximum length of sequence; maxLen is enforced before trimming and truncation; default=Inf (no maximum)
      --minLen                      integer. minLen is enforced after trimming and truncation; default=50
      --rmPhiX                      {"T","F"}. remove PhiX from read              
      --minOverlap                  integer. minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
      --maxMismatch                 integer. The maximum mismatches allowed in the overlap region; default=0
      --trimOverhang                {"T","F"}. If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off. "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction primer region. Default="F" (false)

    Other arguments:
      --dadaOpt.XXX                 Set as e.g. --dadaOpt.HOMOPOLYMER_GAP_PENALTY=-1 Options for the dada function, see ?setDadaOpt in R for available options and their defaults
      --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are  pseudo pooling: "pseudo", true: "T", false: "F"
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random  mnemonic.

     Help:
      --help                        Will print out summary above when executing nextflow run h3abionet/TADA                                   

     Example run:
     To run on UCT hex
     1) Start a 'screen' session from the headnode
     2) Start an interactive job using: qsub -I -q UCTlong -l nodes=1:series600:ppn=1 -d `pwd`
     3) A typical command would look something like:

        nextflow run h3abionet/TADA --trimFor 24 --trimRev 25 --reference /specify/relevant/directory/gg_13_8_train_set_97.fa.gz --email katieviljoen@gmail.com -profile uct_hex --reads  '/specify/relevant/directory/*{R1,R2}.fastq' -with-singularity /scratch/DB/bio/singularity-containers/1a32017e5935-2018-05-31- db3a9cebe9fc.img --pool 'pseudo'


## Prerequisites

Nextflow (>=20.11.0), dada2 (>= 1.8), R (>= 3.2.0), Rcpp (>= 0.11.2), methods (>= 3.2.0), DECIPHER, phangorn, biomformat
Note: if you are working on UCT hex you can simply use the singularity image specified in the uct_hex profile (no need to install these R packages)

## Documentation
The h3abionet/TADA pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)

## Built With

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/what-docker)
* [Singularity](https://sylabs.io/docs/)


## Credits

The initial implementation of the DADA2 pipeline as a Nextflow workflow (https://github.com/HPCBio/dada2-Nextflow) was done by Chris Fields from the high performance computational biology unit at the University of Illinois (http://www.hpcbio.illinois.edu). Please remember to cite the authors of [DADA2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/) when using this pipeline. Further development to the Nextflow workflow and containerisation in Docker and Singularity for implementation on UCT's HPC was done by Dr Katie Lennard and Gerrit Botha, with inspiration and code snippets from Phil Ewels http://nf-co.re/

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
