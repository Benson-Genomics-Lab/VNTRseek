# Overview

VNTRseek is a computational pipeline for the detection of VNTRs. Given an input reference set of TRs and a set of next-generation sequencing reads, such as those produced by Illumina, VNTRseek produces a database and a VCF file with VNTR calls.

Currently, input reads in FASTA and FASTQ format are supported. Files must be compressed either using plain gzip (.gz) or gzip then tar (.tar.gz). Preliminary BAM support is available beginning with version 1.09.

VNTRseek is also capable of running on clusters, though only the Grid Engine queuing system has been tested. Support for other platforms as pull requests is welcome.

# Documentation

The full documentation is located on the [wiki](https://github.com/yzhernand/VNTRseek/wiki) on our [GitHub page](https://github.com/yzhernand/VNTRseek).

# Download

The latest version of VNTRseek is 1.09.2. You can download it from the [releases tab](https://github.com/yzhernand/VNTRseek/releases) on our GitHub repo, or from our [download page](http://orca.bu.edu/vntrseek/download.php).

The master branch is currently tracking the development of the unstable version 1.10.