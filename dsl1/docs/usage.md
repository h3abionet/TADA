# uct-cbio/16S-rDNA-dada2-pipeline Usage

## General Nextflow info
Nextflow handles job submissions on different environments (PBS in our case), and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
For University of Cape Town users who will be running Nextflow on UCT's HPC (hex), you also need to include the following lines in your `~/.bashrc`:

```bash
export JAVA_HOME=/opt/exp_soft/java/jdk1.8.0_31/
JAVA_CMD=/opt/exp_soft/java/jdk1.8.0_31/bin/java
export PATH=$PATH:/opt/exp_soft/cbio/nextflow
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run uct-cbio/16S-rDNA-dada2-pipeline --reads '*_R{1,2}.fastq.gz' --trimFor 24 --trimRev 25 --reference 'gg_13_8_train_set_97.fa.gz' -profile uct_hex
```

This will launch the pipeline with the `uct_hex` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
outdir         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
For more information on dada2 usage and how to import the output from this pipeline into R for downstream analyses using the phyloseq package see https://benjjneb.github.io/dada2/tutorial.html

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull uct-cbio/16S-rDNA-dada2-pipeline
```
Note that `nextflow pull` by default stores pipelines in $HOME/.nextflow/assets To clone this pipeline to a specific directory, use `nextflow clone ct-cbio/16S-rDNA-dada2-pipeline target-dir ` where target-dir is the target directory.

## Main (required) arguments

### `-profile`
Use this parameter to choose a configuration profile. Each profile is designed for a different compute environment - follow the links below to see instructions for running on that system. Available profiles are:

* `standard`
    * The default profile, used if `-profile` is not specified at all. Runs locally and expects all software to be installed and available on the `PATH`.
    * This profile is mainly designed to be used as a starting point for other configurations and is inherited by most of the other profiles.
* `uct_hex`
    * Designed to run on UCT's high-performance cluster (hex).
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--trimFor` and `--trimRev`

Set length of R1 (--trimFor) and R2 (--trimRev) that needs to be trimmed (set 0 if no trimming is needed)

```bash
--trimFor 24 --trimRev 25
```

## Reference Genomes

Specify path to taxonomic database to be used for annotation step

### `--reference` 

```bash
--reference ./gg_13_8_train_set_97.fa.gz
```
See https://benjjneb.github.io/dada2/training.html for dada2-maintained reference files for download

## Full list of arguments (example)
Reads          : data/*{1,2}.fastq.gz
trimFor        : 0
trimRev        : 0
truncFor       : 0
truncRev       : 0
truncQ         : 2 (default)
maxEEFor       : 2
maxEERev       : 2
maxN           : 0 (default)
maxLen         : Inf (default)
minLen         : 20 (default)
rmPhiX         : FALSE (default)
minOverlap     : 20 (default is 12)
maxMismatch    : 0 (default)
trimOverhang   : FALSE
species        : false
pool           : false
Reference      : ./gg_13_8_train_set_97.fa.gz
Max Memory     : 256 GB
Max CPUs       : 64
Max Time       : 24d 20h 31m 24s
Output dir     : ./2018-06-04-dada2
Working dir    : /researchdata/katie/dada2-test/work
Container      : /DB/singularity-containers/1a32017e5935-2018-05-31-db3a9cebe9fc.img
Pipeline Release: dev

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests
Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files in [`conf`](../conf) for examples.

## General line parameters
### `--outdir`
The output directory where the results will be saved.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to speicfy this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the summary HTML / e-mail.

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override defaults. For example, if you don't want FastQC errors to be ignored, you can specify a config file using `-c` that contains the following:

```groovy
process.$fastqc.errorStrategy = 'terminate'
```

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--clusterOptions`
Submit arbitrary cluster scheduler options (not available for all config profiles). For instance, you could use `--clusterOptions '-p devcore'` to run on the development node (though won't work with default process time requests).


---

[![UCT CBIO](/conf/cbio_logo.png)](http://www.cbio.uct.ac.za/)

---
