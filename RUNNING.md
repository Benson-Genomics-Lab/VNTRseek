# How to

1. [Run VNTRseek](#Run-VNTRseek)
   1. [v2.0.3](#v2.0.3)
      - [the Long Way](#the-Long-Way)
      - [the Short Way](#the-Short-Way)
   1. [v2.0.1](#v2.0.1)
   1. [v2.0.0](#v2.0.0)
1. [Specify Options](#Specify-Options)
      - [on the Command Line](#on-the-Command-Line)
      - [in a Configuration File](#in-a-Configuration-File)
      - [for NPROCESSES](#for-NPROCESSES)
1. [Use qsub](#Use-qsub)
1. [Use Cram Files](#Use-Cram-Files)

## Run VNTRseek

You'll need an environment with samtools, Perl (>= 5.24) and SQLite (>= 3.37).

VNTRseek takes as input fasta, fastq, sam, bam, and cram formatted files.
Files in the folder provided to the `--INPUT_DIR` option are processed.
A database is produced, along with several other output files.
Older versions leave many intermediate files; versions >= 2.0.2 allow
for the intermediate files to be deleted on a successful run.

Different versions have different option handling, per the latest development.
Make sure you're reading the section for the relevant version.

VNTRseek must be run with multiple processes (altered via the `--NPROCESSES` option).
This is different from _processors_; there is internal adjustment necessary
for parallelized tasks, see the options sections [for NPROCESSES]($for-NPROCESSES)
for more.

### Input_Dir

Only one format is processed for a run.
The input reader searches `--INPUT_DIR` for each format and uses the first one found.
The current search order is `fasta, fastq, sam, bam, cram`.
Fastas take both `fa` and `fasta` extension;
similarly fastqs include both `fq` and `fastq`.

For example, `--INPUT_DIR` contains three `.fa`s, five `.fasta`s,
two `.fastq`s, and four `.cram`s; the eight fastas would be processed.

Further, if the first (alphabetically) processed file is compressed, it
will attempt to decompress all processed files accordingly.
Supported compressions are `targz gz tarbzip bzip tarxz xz`.
Mixed compressions will result in an error.


## v2.0.3

The minimum options required for running are

- the `--RUN_NAME`,
  - a name used for output files (_not_ a file prefix)
- the `--REFERENCE`,
  - file prefix to reference files for input sequences
- and `--INPUT_DIR` options.
  - directory containing fasta or bam files
  - all such in the folder will be analyzed

Paths can be relative. Be sure to consider [more options](#Specify-Options),
particularly `OUTPUT_DIR` (`.` by default) and `NPROCESSES`.

A database `$RUN_NAME.db` and other output files (including a terse config)
will be produced in
```
$OUTPUT_DIR/vntr_$RUN_NAME/
```

### the Long Way
```sh
path/to/vntrseek.pl --HELP
path/to/vntrseek.pl 0 20 --RUN_NAME "example" --REFERENCE "exGenome/prefix" --INPUT_DIR "exFastas/"
path/to/vntrseek.pl 100 --RUN_NAME "example" # this clears any error codes
```

Above are the typical commands to run from a terminal.
By default it will attempt to use 2 processes (the parent and 1 child for parallelized tasks).
Be mindful of [other options](#Specify-Options) as well.
(Running the help command is not necessary, but recommended).

A start and end step (inclusive) must be provided to run the pipeline itself,
beginning at 0 and ending with 20 for the full pipeline.

VNTRseek stores the last error it encountered and will not resume for a particular
run unless the error is cleared via a step instruction of 100. This diverts from the 
pipeline and exits, so no other steps can be specified.
It needs all information necessary to locate the run database; this is `--RUN_NAME` and,
if you specified one, the `--OUTPUT_DIR`.


### the Short Way

```sh
path/to/vntrseek.pl 0 20 --CONFIG "exFolder/config.cnf"
path/to/vntrseek.pl 100 --CONFIG "exFolder/config.cnf"
```

The `--CONFIG` option can take a file in INI format with any desired running options.
All minimum options must still be provided between the config and command line.
See how to [specify options in a configuration file](#in-a-Configuration-File) for an example.
The `default.vs.cnf` file also lists all running options with verbose descriptions.

You can use the `--GEN_CONFIG` option alone with a path (file or directory).

```sh
.../vntrseek.pl --GEN_CONFIG "exFolder/config.cnf"
```

This generates a default config at the specified path.
Given a directory (that exists), it will create `directory/example.vs.cnf`.
Given a file name (that doesn't exist), it will write the default config in that file.
It will not create any part of the path name that does not exist.

Other options specified at the command line _will override anything_ in the config file.
This can be used as a way of consolidating common options in the config and varying only
what's needed at the command line: e.g. `--INPUT_DIR` and `--RUN_NAME`.
We do not recommend altering the underlying `default.vs.cnf` unless you're
certain you know everything it affects.


## Specify Options

Options can be specified on the command line or from a configuration file.
An option can be specified in both places, in which case
the command line value will be used.

The config path is specified with the `--CONFIG` option,
and can contain _all_ information for running. Paths can be relative, and `--GEN_CONFIG`
produces a default config at a specified path.


### on the Command Line

```
# From --HELP message

Usage: .../vntrseek.pl <start step> <end step> [options]

        The first step is 0, last step is 20.
        At least --RUN_NAME, --INPUT_DIR, and --REFERENCE must be provided,
          or a valid --CONFIG file specifying them.
          Use --GEN_CONFIG to generate a default file.

        A config can be provided with command line arguments,
          in which case the command line values will take precedence.


OPTIONS:

        --HELP                        prints this help message
        --RUN_NAME                    prefix for database name (such as the name of your analysis)
        --INPUT_DIR                   input data directory (BAM or plain or gzipped fasta/fastq files)
        --OUTPUT_DIR                  output directory (must be writable and executable!)
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

        --NPROCESSES                  number of processors to use on your system (default, minimum 2)
        --CONFIG                      path to file with any of the above options. Command line option only.
        --GEN_CONFIG                  set with a path to generate a clean config. Command line option only.
        --CLEANUP                     clear intermediate files if run successful. Command line option only.
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

Any number of options can be passed. All options on the command line will override
any in the config passed to `--CONFIG`. The complete, final set of options will be recorded in a terse
config file in the ouput folder.


### in a Configuration File

The config path is passed to `--CONFIG <config_path>`,
and can contain all options for running.
A config can be generated using the `--GEN_CONFIG <config_path>` option
and a separate, complete configuration of every
run will be generated in the `$OUTPUT_DIR/vntr_$RUN_NAME/` folder.

Example Partial Config:
```sh
REFERENCE=path/to/reference/base_name
INPUT_DIR=path/to/folder/with/fastas_or_bams/
OUTPUT_DIR=path/for/output/

PLOIDY=2
IS_PAIRED_READS=1
READ_LENGTH=150
MIN_FLANK_REQUIRED=10
MAX_FLANK_CONSIDERED=50

NPROCESSES=8
```

See `defaults.vs.cnf` for universal underlying defaults and verbose descriptions.
Modification of this file is not recommended.

### for NPROCESSORS

Perl's `fork` facility does not respect environment variables limiting resources.
So 9 processes will try to use 9 processors, even if the environment only allows 8.
A parent and 8 children will use 9 processes.
For parallelized tasks, the parent waits on children, freeing its processor
while it waits. There will be small windows in which 9 processors are used,
but largely only 8 will be active.


Specifying 8 processes
will guarantee only 8 processors are ever used; 7 children will be
spawned during paralellized tasks.

This difference is important when operating on clusters where jobs are managed by a
grid engine (this is you if you use `qsub`).
Jobs that exceed their resource limits are typically reaped.
Thus if you specify a parallel environment with 8 cores (processors),
and if VNTRseek uses 8 children,
there could be a window in which 9 processors are active,
and the job could be reaped.

For non-managed cases, specifying more processes than the number of cores is not an issue,
so if you'd like to run 8 children on an 8 core machine and there's no overhead threatening
to reap the process group, specify 9 for `--NPROCESSES` in newer versions.


## Use qsub

When using a computer cluster with a grid engine, you'll likely want to use `qsub`
with a script containing the options and environment for VNTRseek.
Below is an example of such a script. A script is submitted to qsub as `qsub myscript.sh`.

```sh
#!/bin/bash -l

#$ -P mylab
#$ -N vntrseek_xmpl
#$ -o /somewhere/qsub.xmpl.out
#$ -e /somewhere/qsub.xmpl.err
#$ -M your.email@domain
#$ -m bea
#$ -l buyin
#$ -l h_rt=24:00:00
#$ -pe omp 8

# your system may have some form of environment manager for
#  loading and unloading tools
module load samtools perl/5.28.1 sqlite3

# specifying 100 clears any error
/path/to/vntrseek.pl 100 --CONFIG xmpl.cnf
# otherwise, specify start and stop steps. 0-20 is the whole pipeline
/path/to/vntrseek.pl 0 20 --CONFIG xmpl.cnf
```

The Sun Grid Engine has a variety of options for running
that can be specified at the command line or in the script on lines starting with `#$`.
Here are a few relevant ones.

```
-P        Name of a project, may be used for determining
            available resources (CPUs, GPUs, RAM)
-N        Name of the job, appears when looking at the queue.

-o        Name for a file that will receive standard output text.
-e        Name for a file that will receive standard error text.
-j        Flag indicating to join the output and error text.
           Default behavior is to not join.
           Adding this option joins the text.
           A y/n value may be specified to be explicit.

-M        Email address alerted on conditions listed by -m option.
-m        Conditions on which to send email to address listed at -M option.
            Here, when the job (b)egins, (e)nds, or is (a)borted.

-l        Variables indicating desired resources.
            buyin means use resources reserved for the project
              specified by the -P option.
            h_rt=HH:MM:SS sets a hard limit on the wall-clock time
              for the job to run.
-pe       Terms indicating a parallel environment to use.
            These will be platform dependent; above, omp is a simple
            multiprocessing environment and we request 8 processors.
```

## Use Cram Files

The samtools `view` command reads cram files based
on the metadata within the cram header.
Since VNTRseek uses this utility,
cram files are intrinsically supported
so long as the environment variables `REF_PATH` and `REF_CACHE` are correct.

See [the samtools doc](https://www.htslib.org/workflow/cram.html#the-ref_path-and-ref_cache)
for information on the `REF_PATH` and `REF_CACHE` environment variables.
If you downloaded a cram without its reference files, be aware that
samtools may attempt an internet connection to retrieve them.
This process may result in a large file cache being created
and significant slow down based on the network.