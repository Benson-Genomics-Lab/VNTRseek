# Overview

VNTRseek is a computational pipeline for the detection of VNTRs.
Given an input reference set of TRs and a set of next-generation sequencing reads,
such as those produced by Illumina, VNTRseek produces a database and a VCF file with VNTR calls.

Currently, input reads in FASTA, FASTQ, SAM, BAM, and CRAM format are supported.
Files can be compressed using plain gzip (.gz), bzip2 (.bz, .bz2), or xz (.xz) format,
optionally preceeded by tar (.tgz, .tar.gz, .tbz, .tbz2, .tar.bz , .tar.bz2, .txz, .tar.xz)

VNTRseek readily runs on Sun Grid Engine clusters. Support for other platforms as pull requests is welcome.

# Documentation

See the [running instructions](RUNNING.md) for information on the latest running parameters.
An [output description](RESULTS.md) is included.

# Download

The latest version of VNTRseek is 2.0.3. You can download it from the [Github repo](https://github.com/Benson-Genomics-Lab/VNTRseek/releases).
An active (likely unstable) [branch for v2 series development](https://github.com/Benson-Genomics-Lab/VNTRseek/tree/v2.0_dev) is available for pull contributions.

# Edlib Attribution

A slimmed down version of [Edlib](https://github.com/Martinsos/edlib) is included.
It contains the original code for calculations, all credit to original author.
Test scripts and extended features are removed, as well as original build files.
A wrapper script has been added for obtaining the NICE alignment
strings from simple command line input.