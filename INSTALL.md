# Important Notes

- Only CentOS 6 and Ubuntu 12.10 and up have been tested.
- Only a Linux 64-bit version of our software is available at this time.
- Currently, the default install path is valid only on UNIX-like platforms.
- Version 1.10 and later has switched to using SQLite for database management.


# Requirements

## Hardware

VNTRseek was designed around multiprocessing on a large computer.
Minimum specs suggested are 8 cores, 32GB of RAM, and 1TB of file space.
This will run a sample of 700M human FASTA reads in roughly 30 hours.

## Software

Latest versions preferred.

- GCC (>= 4.1.2)
- Perl (>= 5.24.0)
  - DBI and DBD::SQLite modules
- SQLite (>= 3.37.0)
- TRF (>= 4.09)
- samtools (>= 1.8)
- GLIBC (>= 2.14)
- CMake (>= 3.2)

TRF is downloaded during installation. If the download fails, you can download it manually from
[github](https://github.com/Benson-Genomics-Lab/TRF/releases/latest/download)
and save it as `trf409-ngs.linux.exe` in the build directory (see [Installation](#installation)
below).

# Installation

To build and install to the default directory, download the `VNTRseek-N.N.N.tar.gz` file from the latest release and run the following commands:

```sh
tar xzvf VNTRseek-N.N.N.tar.gz
cd VNTRseek-N.N.N
mkdir build
cd build
cmake ..
make install # or sudo make install, if needed
```

By default, this will install the pipeline to `/usr/local/vntrseekN.N.N` (e.g.,
`/usr/local/vntrseek2.0.3`) and create a symbolic link at `/usr/local/bin/vntrseek` 
that points to the main program at `/usr/local/vntrseekN.N.N/vntrseek.pl`.

If you would like to choose a different installation prefix,
add the `-DCMAKE_INSTALL_PREFIX` option to the `cmake` call, e.g.:

```sh
# use this for an absolute path to vntrseekN.N.N
cmake -DCMAKE_INSTALL_PREFIX=<absolute_prefix> ..

# use this to create a subdirectory under your home directory for vntrseekN.N.N
cmake -DCMAKE_INSTALL_PREFIX=${HOME} ..
```

If you install to a non-standard location, you may need to add the `/bin` with the link
(e.g., `<absolute_prefix>/bin` or `${HOME}/bin`) to your `PATH` variable.

You can change your default compiler using the `-DCMAKE_C_COMPILER` option.
For example, to use clang:

```sh
cmake -DCMAKE_C_COMPILER=clang ..
```

**If you installed this pipeline as root, and are creating an INDIST
file** you may need to run it as root unless you give your user
permission to write to the installation directory.


# Uninstalling

On UNIX, run:

```sh
xargs rm < install_manifest.txt # or sudo xargs rm < install_manifest.txt
```

from the build directory you created above. The directory will
remain, however, so you will not lose any reference files.
