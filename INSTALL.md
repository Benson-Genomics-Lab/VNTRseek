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

- samtools
- Perl (>= 5.24.0)
  - DBI and DBD::SQLite modules
- SQLite (>= 3.37.0)

TRF is also required, but is downloaded during installation. If the download fails,
you can download it manually from
[the website](http://tandem.bu.edu/trf/trf409.linux64.download.html)
and save it as `trf409-ngs.linux.exe` in the build directory (see [Installation](#installation)
below).

Building requires cmake (http://www.cmake.org/), minimum version 3.0

Additionally you will need GCC version 4.1.2 on Mac/Linux/CYGWIN or
a compatible compiler.


# Installation

To build and install to the default directory, run the following commands:

```sh
tar xzvf vntrseekN.NN.tar.gz
cd vntrseekN.NNsrc
mkdir build
cd build
cmake ..
make install # or sudo make install, if needed
```

By default, this will install the pipeline to `/usr/local/vntrseekN.NN` (eg,
`/usr/local/vntrseek2.0.0`).

If you would like to choose a different installation prefix,
add the `-DCMAEK_INSTALL_PREFIX` option to the `cmake` call, e.g.:

```sh
cmake -DCMAKE_INSTALL_PREFIX=<absolute_prefix> ..

# install to home at ${HOME}/vntrseekN.NN
cmake -DCMAKE_INSTALL_PREFIX=${HOME} ..
```

If you install to a non-standard location, you may need to add
`<absolute_prefix>/bin` to your `PATH` variable (e.g., if your prefix
was `/opt`, you will need to have `/opt/bin` in your `PATH`).

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
