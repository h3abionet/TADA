<p>
<img align="left" src="./assets/cbio_logo.png" width="200" hspace="50"/>
<img align="left" src="./assets/HPCBio.png" width="200" hspace="50"/>
</br></br></br>
</p>

# TADA - Targeted Amplicon Diversity Analysis using DADA2, implemented in Nextflow

A dada2-based workflow using the Nextflow workflow manager for Targeted Amplicon Diversity Analysis.

## Badges
| fair-software.nl recommendations                        |                             |
| ------------------------------------------------------- | --------------------------- |
|(1/5) code repository                                    |[![GitHub Repo Status](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/h3abionet/TADA)|
|(2/5) license                                            |[![GitHub License Status](https://img.shields.io/github/license/h3abionet/TADA)](https://github.com/h3abionet/TADA)|
|(3/5) community registry                                 | [bio.tools Registry](https://bio.tools/tada-amplicon) |
|(4/5) citation                                           |[![Zenodo Status](https://zenodo.org/badge/DOI/10.5281/zenodo.4208836.svg)](https://doi.org/10.5281/zenodo.42088362)|
|(4/5) checklist                                          | |
|overall                                                  |[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)|
**GitHub Actions**
|Docker build                                             | [![GitHub Docker Status](https://github.com/h3abionet/TADA/actions/workflows/docker.yml/badge.svg)](https://github.com/h3abionet/TADA/actions?query=workflow%3A%22Docker%22)|
|Continuous integration                                   | [![GitHub CI Status](https://github.com/h3abionet/TADA/actions/workflows/ci.yml/badge.svg)](https://github.com/h3abionet/TADA/actions?query=workflow%3A%22CI%22)|

## Basic usage:

The latest help menu can be accessed using `nextflow run h3abionet/TADA --help`.  

```
  Usage:

  This pipeline can be run specifying parameters in a config file or with command line flags.
  The typical example for running the pipeline with command line flags is as follows:
  
    nextflow run h3abionet/TADA --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 \
      --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex

  The typical command for running the pipeline with your own config (instead of command line flags) is as follows:
  
    nextflow run h3abionet/TADA -c dada2_user_input.config -profile uct_hex
  
  where 'dada2_user_input.config' is the configuration file (see example 'dada2_user_input.config')
  
  NB: '-profile uct_hex' still needs to be specified from the command line

  Parameters
  ----------

  Mandatory arguments:
    -profile                      Hardware config to use. Currently profile available for UCT's HPC 'uct_hex' and UIUC's 'uiuc_singularity' - create your own if necessary
                                  NB -profile should always be specified on the command line, not in the config file      
    
  Input (mandatory): Additionally, only one of the following must be specified:
    --reads                       Path to FASTQ read input data.  If the data are single-end, set '--single-end' to true.
    --input                       Path to a sample sheet (CSV); sample sheet columns must have a headers with 'id,fastq_1,fastq_2'.  
    --seqTables                   Path to input R/dada2 sequence tables. Only sequence tables with the original ASV sequences as the identifier are supported

  Output location:
    --outdir                      The output directory where the results will be saved

  Read preparation parameters:
    --trimFor                     integer. Headcrop of read1 (set 0 if no trimming is needed)
    --trimRev                     integer. Headcrop of read2 (set 0 if no trimming is needed)
    --truncFor                    integer. Truncate read1 here (i.e. if you want to trim 10bp off the end of a 250bp R1, truncFor should be set to 240). Enforced before trimFor/trimRev
    --truncRev                    integer. Truncate read2 here ((i.e. if you want to trim 10bp off the end of a 250bp R2, truncRev should be set to 240). Enforced before trimFor/trimRev
    --maxEEFor                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
    --maxEERev                    integer. After truncation, R1 reads with higher than maxEE "expected errors" will be discarded. EE = sum(10^(-Q/10)), default=2
    --truncQ                      integer. Truncate reads at the first instance of a quality score less than or equal to truncQ; default=2
    --maxN                        integer. Discard reads with more than maxN number of Ns in read; default=0
    --maxLen                      integer. Maximum length of trimmed sequence; maxLen is enforced before trimming and truncation; default=Inf (no maximum)
    --minLen                      integer. Minimum length enforced after trimming and truncation; default=50
    --rmPhiX                      {"T","F"}. remove PhiX from read

    In addition due to modifications needed for variable-length sequences (ITS), the following are also supported.  Note if these are set,
    one should leave '--trimFor/--trimRev' set to 0.

    --fwdprimer                   Provided when sequence-specific trimming is required (e.g. ITS sequences using cutadapt).  Experimental
    --revprimer                   Provided when sequence-specific trimming is required (e.g. ITS sequences using cutadapt).  Experimental

  Read merging:
    --minOverlap                  integer. minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
    --maxMismatch                 integer. The maximum mismatches allowed in the overlap region; default=0
    --trimOverhang                {"T","F"}. If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off.
                                  "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen when reads are longer than the amplicon and read into the other-direction                                               primer region. Default="F" (false)

  Error models:
    --qualityBinning              Binned quality correction (e.g. NovaSeq/NextSeq).  default: false
    --errorModel                  NYI. Error model to use (one of 'illumina', 'illumina-binned', 'pacbio-ccs', 'custom'). This will replace
                                  '--qualityBinning'

  Denoising using dada:
    --dadaOpt.XXX                 Set as e.g. --dadaOpt.HOMOPOLYMER_GAP_PENALTY=-1 Global defaults for the dada function, see ?setDadaOpt in R for available options and their defaults
    --pool                        Should sample pooling be used to aid identification of low-abundance ASVs? Options are
                                  pseudo pooling: "pseudo", true: "T", false: "F"

  Merging arguments (optional):
    --minOverlap                  The minimum length of the overlap required for merging R1 and R2; default=20 (dada2 package default=12)
    --maxMismatch                 The maximum mismatches allowed in the overlap region; default=0.
    --trimOverhang                If "T" (true), "overhangs" in the alignment between R1 and R2 are trimmed off. "Overhangs" are when R2 extends past the start of R1, and vice-versa, as can happen
                                  when reads are longer than the amplicon and read into the other-direction primer region. Default="F" (false)
    --minMergedLen                Minimum length of fragment *after* merging; default = 0 (no minimum)
    --maxMergedLen                Maximum length of fragment *after* merging; default = 0 (no maximum)

  ASV identifiers:
    --idType                      The ASV IDs are renamed to simplify downstream analysis, in particular with downstream tools.  The
                                  default is "md5" which will run MD5 on the sequence and generate a QIIME2-like unique hash.  Alternatively, 
                                  this can be set to "ASV", which simply renames the sequences in sequencial order.  

  Taxonomic arguments.  If unset, taxonomic assignment is skipped
    --taxassignment               Taxonomic assignment method.  default = 'rdp'
    --reference                   Path to taxonomic database to be used for annotation (e.g. gg_13_8_train_set_97.fa.gz). default = false
    --species                     Specify path to fasta file. See dada2 addSpecies() for more detail. default = false
    --minBoot                     Minimum bootstrap value.  default = 50
    --taxLevels                   Listing of taxonomic levels for 'assignTaxonomy'. Experimental.

  Chimera detection:
    --skipChimeraDetection        Skip chimera detection/removal; default = false
    --removeBimeraDenovoOpts      Additional removeBimeraDenovo options; default = ''

  ASV multiple sequence alignment:
    --skipAlignment               Skip alignment step; note this also skips ML phylogenetic analysis. default = false
    --aligner                     Aligner to use, options are 'DECIPHER' or 'infernal'. default = 'DECIPHER'
    --infernalCM                  Covariance model (Rfam-compliant) to use.  default = false.

  Phylogenetic analysis:
    --runTree                     Tool for ML phylogenetic analysis.  Options are 'phangorn' and 'fasttree'. default = 'phangorn'

  Additional output:
    --toBIOM                      Generate a BIOM v1 compliant output. default = true
    --toQIIME2                    Generate QZA artifacts for all data for use in QIIME2. default = false

  Sample names:
    --sampleRegex                 Modify sample names based on a regular expression. default = false.  Note this option
                                  is deprecated in favor of using a sample sheet.

  Additional options:
    --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run
                                  sent to you when the workflow exits
    -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

  Help:
    --help                        Will print out summary above when executing nextflow run uct-cbio/16S-rDNA-dada2-pipeline
```

## Prerequisites

Nextflow (>=20.11.0) with either Singularity (>3.4.1) or Docker (>20.10.1, though we recommend the latest based on security updates).  We don't directly support non-containerized options (locally installed tools, bioconda) though these could work
based on configuration.

## Documentation

The h3abionet/TADA pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)

## Built With

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://www.docker.com/what-docker)
* [Singularity](https://sylabs.io/docs/)

## Credits

The initial implementation of the DADA2 pipeline as a Nextflow workflow (https://github.com/HPCBio/16S-rDNA-dada2-pipeline) was done by Chris Fields from the High Performance Computating in Biology group at the University of Illinois (http://www.hpcbio.illinois.edu). Please remember to cite the authors of [DADA2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927377/) when using this pipeline. Further development to the Nextflow workflow and containerisation in Docker and Singularity for implementation on UCT's HPC was done by Dr Katie Lennard and Gerrit Botha, with inspiration and code snippets from Phil Ewels http://nf-co.re/

## Contributors

The following have contributed to the development, testing, and deployment of this workflow. For the most up-to-date listing see the [Contributors](https://github.com/h3abionet/TADA/graphs/contributors) link.

* [Katie Lennard](https://github.com/kviljoen)
* [Gerrit Botha](https://github.com/grbot)
* [Chris Fields](https://github.com/cjfields)
* [Jessica Holmes](https://github.com/jrkirk61)
* [Gloria Rendon](https://github.com/grendon)
* [Lindsay Clark](https://github.com/lvclark)
* [Wojtek Bażant](https://github.com/wbazant)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
