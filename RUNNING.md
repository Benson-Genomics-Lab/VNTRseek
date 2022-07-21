# How to

1. [Run VNTRseek v2.0.0](#Run-VNTRseek-v2.0.0)
    - [the Long Way](#the-Long-Way)
    - [the Short Way](#the-Short-Way)
1. [Run VNTRseek v2.0.1](#Run-VNTRseek-v2.0.1)
1. [Specify Options](#Specify-Options)
    - [on the Command Line](#on-the-Command-Line)
    - [in a Configuration File](#in-a-Configuration-File)
1. [Use qsub](#Use-qsub)

## Run VNTRseek v2.0.0

### the Long Way
```sh
# at terminal
path/to/vntrseek.pl 100 # this clears any error codes
path/to/vntrseek.pl 0 19 --DBSUFFIX "example" --REFERENCE "exGenome" --INPUT_DIR "exFastas"
  --OUTPUT_ROOT "exOut" --TMPDIR "exTmp" --SERVER "ex.am.ple" --NPROCESSES 1
```

The minimum required for running are the `--DBSUFFIX`, `--REFERENCE`, `--INPUT_DIR`,
and `--OUTPUT_ROOT`options. All paths must be absolute.

The pipeline will automatically look for [a configuration file](#in-a-Configuration-File) with the name `<DBSUFFIX>.vs.cnf`
in the directory where vntrseek was called that can contain any options besides `--DBSUFFIX`.

The `--REFERENCE` flag takes a prefix (no file extension) for the desired reference genome (fasta format).
A few further files are needed, a `.leb36` and `.seq` file as described in REFERENCE_FORMATS.

VNTRseek takes as input fasta and bam formatted files. All such files in the folder
provided to `--INPUT_DIR` are processed.

The `--OUTPUT_ROOT` is a path where you want the output generated. Similarly `--TMPDIR`
is for temporary output that is meant to be deleted after a run. Should a run fail,
these files can be safely deleted; they are _atomic_ to the step that generates them.

The `--SERVER` is for writing absolute html link addresses in the output.

Specify a number greater than 1 for `--NPROCESSES`, a number of
processing cores available in the execution environment.

### the Short Way

At the command line, you'll enter something like

```sh
# at terminal
path/to/vntrseek.pl 100
path/to/vntrseek.pl 0 19 --DBSUFFIX "example"
```

and you will need to [Specify Options in a Configuration File](#in-a-Configuration-File)

## Run VNTRseek v2.0.1

Read the running instructions for version 2.0.0 first.

Now, version 2.0.1 is run very similarly, except it doesn't look for a config
file based on `--DBSUFFIX`. To supply a config file, use the `--CONFIG`
option with the path to the config. A default config can be generated with
`--GEN_CONFIG` (provided a path). A terse config is output with every run at

```$OUTPUT_ROOT/vntr_$DBSUFFIX/$DBSUFFIX.vs.cnf```

with all final options for the sake of reproducibility.

As a consequence of the changes, `--DBSUFFIX` can now be specified in the config.
Command line options still take precedence over config options. The terse config
will always contain the latest settings used for the output and dbsuffix.

Also, paths no longer need to be absolute and a default path for
`OUTPUT_ROOT` has been added ("output").

## Specify Options

Options can be specified on the command line or from a configuration file.
An option can be specified in both places, in which case
the command line value will be used.

The minimum required for running are the `--DBSUFFIX`, `--REFERENCE`, `--INPUT_DIR`,
and `--OUTPUT_ROOT` options for v2.0.0 and all paths must be absolute.

For v2.0.1, paths can be relative and a default for `OUTPUT_ROOT` has been added.

### on the Command Line

The description here applies to v2.0.1, though most of the options are also
available for v2.0.0, see `--HELP` option.

```
# From --HELP message

Usage: .../vntrseek2.0.1/vntrseek.pl <start step> <end step> [options]

        The first step is 0, last step is 19.
        At least --DBSUFFIX, --INPUT_DIR, and --REFERENCE must be provided,
          or a valid --CONFIG file specifying them.
          Use --GEN_CONFIG to generate a default file.

        A config can be provided with command line arguments,
          in which case the command line values will take precedence.


OPTIONS:

        --HELP                        prints this help message
        --DBSUFFIX                    prefix for database name (such as the name of your analysis)
        --INPUT_DIR                   input data directory (BAM or plain or gzipped fasta/fastq files)
        --OUTPUT_ROOT                 output directory (must be writable and executable!)
        --TMPDIR                      temp (scratch) directory (must be writable!)
        --SERVER                      server name, used for generating html links

        --REFERENCE                   base name of reference files (default set in global config file)
        --REDO_REFDB                  force reinitialization of the reference set database, 0/1 (default 0)
        --REFERENCE_INDIST_PRODUCE    generate a file of indistinguishable references, 0/1 (default 0)

        --PLOIDY                      sample's ploidy (default 2)
        --READ_LENGTH                 length of fasta reads (default 150)
        --IS_PAIRED_READS             data is paired reads, 0/1 (default 1)
        --STRIP_454_KEYTAGS           for 454 platform, strip leading 'TCAG', 0/1 (default 0)
        --KEEPPCRDUPS                 whether to find and remove PCR duplicates. (default: 1, duplicates are kept)

        --MIN_FLANK_REQUIRED          minimum required flank on both sides for a read TR to be considered (default 10)
        --MAX_FLANK_CONSIDERED        maximum flank length used in flank alignments, set high to use full flank (default 50)
        --MIN_SUPPORT_REQUIRED        minimum number of mapped reads which agree on copy number to call an allele (default 2)

        --NPROCESSES                  number of processors to use on your system (default 2)
        --CONFIG                      path to file with any of the above options. Command line option only.
        --GEN_CONFIG                  set with a path to generate a clean config. Command line option only.
        --CLEAN                       force reinitialization of the run database. Command line option only.
        --STATS                       print out a simple table of run statistics. Command line option only.


ADDITIONAL USAGE:

        .../vntrseek.pl 100           clear error
        .../vntrseek.pl 100 N         clear error and set NextRunStep to N (0-19)
                                        This is only when running on a cluster using the
                                        advanced cluster script that checks for NextRunStep

        .../vntrseek.pl --REFERENCE </path/to/reference/basename> [--REDO_REFDB]
          Initialize a reference set. This only needs to be done once for a reference set.
          Supply --REDO_REFDB to force a recreation of the database.
```

### in a Configuration File

In version 2.0.0, if `--DBSUFFIX <run_name>` is provided, all other options can be kept
in `<run_name>.vs.cnf` in the same place where vntrseek gets called.

In version 2.0.1, the config path is passed to `--CONFIG <config_path>`, and can contain
all options for running, minimum `DBSUFFIX`, `INPUT_DIR`, and `REFERENCE`. A config
can be generated using the `--GEN_CONFIG` option and a complete configuration of every
run will be generated in the output folder.

Example Partial Config:
```sh
# In configuration file
REFERENCE=path/to/reference/folder
INPUT_DIR=path/to/fastas/or/bams
OUTPUT_ROOT=path/for/output
TMPDIR=path/for/temp/files

PLOIDY=2
IS_PAIRED_READS=1
READ_LENGTH=150
MIN_FLANK_REQUIRED=10
MAX_FLANK_CONSIDERED=50

NPROCESSES=8
```


## Use qsub

In a file you'll want to write something like the following. Fill in the details
to suit your use. To run the script with qsub, run at the command line `qsub myscript.sh`
(though perhaps give the script a more informative name).

```sh
#!/bin/bash

#$ -P vntrseek
#$ -N vntrseek_xmpl
#$ -o /somewhere/qsub.xmpl.out
#$ -e /somewhere/qsub.xmpl.err
#$ -M your.email@domain
#$ -m bea
#$ -l buyin
#$ -l h_rt=24:00:00
#$ -pe omp 8

module load samtools perl/5.28.1

# specifying 100 clears any error
/path/to/vntrseek.pl 100 --DBSUFFIX xmpl
# otherwise, specify start and stop steps. 0-19 is the whole pipeline
/path/to/vntrseek.pl 0 19 --DBSUFFIX xmpl
```

The SGE has a variety options for running that are specified on the lines starting with `#$`.

```
-P        Name of a project, used for determining
            available resources (CPUs, GPUs, RAM)
-N        Name of the job, appears when looking at the queue.

-o        Name for a file that will receive standard output text.
-e        Name for a file that will receive standard error text.
-j        Flag indicating to join the output and error text.
           Default behavior is to not join.
           Adding this option joins the text.
           A y/n value may be specified to be explicit.

-M        Email address alerted on conditions listed by -m options.
-m        Conditions on which to send email to address listed at -M option.
            Here, when the job (b)egins, (e)nds, or is (a)borted.

-l        Variables indicating desired resources.
            buyin means use resources purchased for the project
              specified by the -P option.
            h_rt=HH:MM:SS sets a hard limit on the wall-clock time
              for the job to run.
-pe       Variables indicating a parallel environment to use.
            These will be platform dependent; here omp is a simple
            multiprocessing environment and we request 8 processors.
```

