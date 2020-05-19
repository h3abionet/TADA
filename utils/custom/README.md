# Custom primer helper scripts

These are some simple utility scripts that might help if you need to assign NCBI taxonomy ranks to a set of data. 

The scripts use perl (>5.16) and BioPerl 1.7.7 or higher, with SQLite 3.7.17.  All tools were run on the IGB Biocluster:

http://biocluster.igb.illinois.edu 

## Load taxonomy database

Simple BioPerl script, requires [NCBI taxonomy](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).  Uncompress and untar the taxonomy files into a working directory and run the script from there.

```
perl load-taxonomy.pl
```

This will generate a `taxonomy.sqlite` local database.

## Load mapping database

Data is retrieved from the NCBI taxonomy folder in the mapping files, using the [nucleotide accession file](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)

The mapping database is generated from standard command line using a SQLite loading script (`mapping.sql`):

```
zcat nucl_gb.accession2taxid.gz | sqlite3 --init mapping.sql acc2taxid.sqlite
```

## Add taxonomic ranks to BLAST results

BLAST results were generated using BLAST+ 2.9.0.  Example script:

```
#!/bin/bash
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -p hpcbio
#SBATCH --mem=48000
#SBATCH -J BLASTN

# Strict mode: http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

# Check yer variables (probably redundant with the above but good practice)
# if [ -z ${var+x} ]; then echo "var is unset"; else echo "var is set to '$var'"; fi

# Load Modules
module load BLAST+/2.9.0-IGB-gcc-4.9.4

# BLAST DATABASE
DB=./mitochondria/mitochondrion

# QUERY FASTA
QUERY=asvs.simple.fna

NAME=$( basename $QUERY )

# Run app on file
blastn -db $DB \
     -query $QUERY \
     -evalue 1e-3 -num_threads $SLURM_NPROCS \
     -outfmt 11 \
     -out $NAME.blastn.asn1

blast_formatter -archive $NAME.blastn.asn1 -outfmt "6 qseqid qlen sseqid slen stitle pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore" -max_target_seqs 10 -out $NAME.blastn.txt
blast_formatter -archive $NAME.blastn.asn1 -outfmt 5 -max_target_seqs 10 -out $NAME.blastn.xml
```

Rank information (along with scientific name and common name data) is appended to the BLAST output using the  `blast-taxid-summary.pl` script:

```
perl blast-taxid-summary.pl 
```

