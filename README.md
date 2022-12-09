# Overview

VNTRseek is a computational pipeline for the detection of VNTRs. Given an input reference set of TRs and a set of next-generation sequencing reads, such as those produced by Illumina, VNTRseek produces a database and a VCF file with VNTR calls.

Currently, input reads in FASTA, FASTQ, and BAM format are supported. Files can be compressed either using plain gzip (.gz) or gzip then tar (.tar.gz).

VNTRseek readily runs on Sun Grid Engine clusters. Support for other platforms as pull requests is welcome.

# Documentation

See the [running instructions](RUNNING.md) for information on the latest running parameters.

# Download

The latest version of VNTRseek is 2.0.2. You can download it from the [Github repo](https://github.com/Benson-Genomics-Lab/VNTRseek/releases).
An active (likely unstable) [branch for v2 series development](https://github.com/Benson-Genomics-Lab/VNTRseek/tree/v2.0_dev) is available for pull contributions.