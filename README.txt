PIGEON: a pipeline for GEO
##########################

INSTALL
-------

Pigeon is a Python 3 script that requires the following programs:
    - Bowtie 2
    - RSEM
    - Samtools
    - STAR

These programs do not need to be in the $PATH, as their binary path can be defined in a configuration file.
However, sort (GNU coreutils) and gzip are expected to be in the $PATH, as they are installed on most UNIX-based systems.

Two C programs have to be compiled. Use the make command in the utils/ directory to generate the binaries.

You will also need to generate aligners indexes to align FASTQ files. Please refer to the aligners manual to know how to generate indexes.


CONFIGURATION
-------------

Pigeon needs a configuration file to run. This file contains :
    - paths to binaries (e.g. Bowtie2)
    - Pigeon API connection details
    - working/output directories
    - alignment/processing details

An empty config.ini file can be found in this directory.


RUNNING PIGEON
--------------

Run Pigeon with the pigeon.py script.

Required options:
    - c FILE            path the configuration file
    - d DATATYPES       one or more data types; available data types are ChIP, HiC, RNA;
                        this option defines which data sets are going to be processed

Ohter options:
    - i FILE            JSON file containing the data sets to process; do not use the API
    - a ASSEMBLIES      genome assemblies to process
    --download INT      number of workers that download SRA files (default: 1)
    --extract INT       number of workers that extract FASTQ files from SRA files (default: 1)
    --align INT         number of workers that align FASTQ files (default: 1)
    --merge INT         number of workers that merge BED files (default: 1)
    --analyze INT       number of workers that analyze merged BED files (default: 1)
    -m, --maxmem INT    buffer size for sorting (in Mb)
    --threads-aln       number of threads that each alignment worker is allowed to use (default: 1)
    --threads-alz       number of threads that each analyze worker is allowed to use (default: 1)


WORKERS
-------

There are five types of worker.
1. Download:    use wget to download a remote SRA file.
2. Extract:     use the SRA fastq-dump to extract FASTQ files from an SRA file.
3. Align:       use Bowtie2 or STAR to align sequences against a reference genome.
4. Merge:       merge and sort data sets
5. Analyze:     run, if implemented, analysis programs on data sets.
                For RNA-seq data sets, RSEM is run to calculare the gene expression.
                For ChIP-seq and Hi-C data sets, nothing is run.




